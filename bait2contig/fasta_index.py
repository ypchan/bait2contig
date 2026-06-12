"""Text FASTA index with metadata needed by bait2contig."""

from __future__ import annotations

import gzip
import mmap
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Optional

from .io import FastaRecord, ensure_parent_dir, fasta_id_from_header, infer_circular_from_header, open_text


SCHEMA_VERSION = "1"
FastaProgressCallback = Callable[[int, int, int], None]


@dataclass(frozen=True)
class IndexedFastaRecord:
    """A FASTA index record without loading the full sequence."""

    id: str
    header: str
    length: int
    is_circular: bool
    seq_offset: Optional[int] = None
    line_bases: int = 0
    line_width: int = 0


def default_index_path(source_path: str | Path) -> str:
    """Return the default text index path for a FASTA input."""

    return f"{source_path}.bait2contig.fai"


def is_gzip_path(path: str | Path) -> bool:
    return str(path).endswith(".gz")


def source_signature(path: str | Path) -> dict[str, str]:
    source = Path(path)
    stat = source.stat()
    return {
        "schema_version": SCHEMA_VERSION,
        "source_path": str(source.resolve()),
        "source_size": str(stat.st_size),
        "source_mtime_ns": str(stat.st_mtime_ns),
        "source_is_gzip": "true" if is_gzip_path(source) else "false",
    }


def index_is_current(index_path: str | Path, source_path: str | Path) -> bool:
    """Return whether an existing text index matches the source FASTA signature."""

    path = Path(index_path)
    if not path.exists() or path.stat().st_size <= 0:
        return False
    try:
        metadata, _ = read_index(path, load_records=False)
    except Exception:
        return False
    expected = source_signature(source_path)
    return all(metadata.get(key) == value for key, value in expected.items())


def open_binary_fasta(path: str | Path):
    if is_gzip_path(path):
        return gzip.open(path, "rb")
    return open(path, "rb")


def build_fasta_index(
    source_path: str | Path,
    index_path: str | Path,
    *,
    progress: Optional[FastaProgressCallback] = None,
    threads: int = 1,
) -> int:
    """Build a text FASTA index and return the number of indexed records."""

    source_path = Path(source_path)
    if not is_gzip_path(source_path):
        return build_plain_fasta_index(source_path, index_path, progress=progress, threads=threads)
    return build_streaming_fasta_index(source_path, index_path, progress=progress)


def build_streaming_fasta_index(
    source_path: str | Path,
    index_path: str | Path,
    *,
    progress: Optional[FastaProgressCallback] = None,
) -> int:
    """Build a text FASTA index with sequential streaming for gzip inputs."""

    source_path = Path(source_path)
    index_path = Path(index_path)
    ensure_parent_dir(index_path)
    tmp_path = index_path.with_name(f"{index_path.name}.tmp.{os.getpid()}")
    tmp_path.unlink(missing_ok=True)

    random_access = not is_gzip_path(source_path)
    signature = source_signature(source_path)
    signature["random_access"] = "true" if random_access else "false"

    records_seen = 0
    total_bases = 0
    seen_ids: set[str] = set()
    header: Optional[str] = None
    record_id: Optional[str] = None
    seq_offset: Optional[int] = None
    sequence_length = 0
    line_bases = 0
    line_width = 0

    with open(tmp_path, "wt", encoding="utf-8", newline="") as out_handle:
        for key, value in signature.items():
            out_handle.write(f"#bait2contig_{key}={value}\n")
        out_handle.write("#columns=id\tlength\toffset\tline_bases\tline_width\tis_circular\theader\n")

        def flush() -> None:
            nonlocal header, record_id, seq_offset, sequence_length, line_bases, line_width, records_seen, total_bases
            if header is None or record_id is None:
                return
            if record_id in seen_ids:
                raise ValueError(f"duplicate FASTA ID in contig FASTA: {record_id}")
            seen_ids.add(record_id)
            out_handle.write(
                "\t".join(
                    [
                        record_id,
                        str(sequence_length),
                        str(seq_offset if seq_offset is not None else -1),
                        str(line_bases),
                        str(line_width),
                        "1" if infer_circular_from_header(header) else "0",
                        header,
                    ]
                )
                + "\n"
            )
            records_seen += 1
            total_bases += sequence_length
            if progress is not None:
                progress(0, 1, 0)
            header = None
            record_id = None
            seq_offset = None
            sequence_length = 0
            line_bases = 0
            line_width = 0

        try:
            with open_binary_fasta(source_path) as in_handle:
                while True:
                    raw_line = in_handle.readline()
                    if not raw_line:
                        break
                    if progress is not None:
                        progress(len(raw_line), 0, 0)
                    stripped = raw_line.rstrip(b"\n\r")
                    if not stripped:
                        continue
                    if stripped.startswith(b">"):
                        flush()
                        header = stripped[1:].decode("utf-8").strip()
                        if not header:
                            raise ValueError("empty FASTA header in contig FASTA")
                        record_id = fasta_id_from_header(header)
                        seq_offset = in_handle.tell() if random_access else None
                        sequence_length = 0
                        line_bases = 0
                        line_width = 0
                        continue
                    if header is None:
                        raise ValueError("sequence found before FASTA header in contig FASTA")
                    sequence = stripped.strip()
                    sequence_length += len(sequence)
                    if line_bases == 0:
                        line_bases = len(sequence)
                        line_width = len(raw_line)
                    if progress is not None:
                        progress(0, 0, len(sequence))
            flush()
            out_handle.write(f"#bait2contig_record_count={records_seen}\n")
            out_handle.write(f"#bait2contig_total_bases={total_bases}\n")
        except Exception:
            tmp_path.unlink(missing_ok=True)
            raise

    os.replace(tmp_path, index_path)
    return records_seen


def build_plain_fasta_index(
    source_path: str | Path,
    index_path: str | Path,
    *,
    progress: Optional[FastaProgressCallback] = None,
    threads: int = 1,
) -> int:
    """Build a text FASTA index for plain FASTA using mmap chunk scanning."""

    source_path = Path(source_path)
    index_path = Path(index_path)
    ensure_parent_dir(index_path)
    tmp_path = index_path.with_name(f"{index_path.name}.tmp.{os.getpid()}")
    tmp_path.unlink(missing_ok=True)

    signature = source_signature(source_path)
    signature["random_access"] = "true"
    threads = max(1, int(threads))

    try:
        with open(source_path, "rb") as in_handle, mmap.mmap(in_handle.fileno(), 0, access=mmap.ACCESS_READ) as data:
            header_starts = find_header_starts(data, threads=threads, progress=progress)
            if not header_starts:
                raise ValueError("contig FASTA has no FASTA records")
            records_seen = 0
            total_bases = 0
            seen_ids: set[str] = set()

            with open(tmp_path, "wt", encoding="utf-8", newline="") as out_handle:
                for key, value in signature.items():
                    out_handle.write(f"#bait2contig_{key}={value}\n")
                out_handle.write("#columns=id\tlength\toffset\tline_bases\tline_width\tis_circular\theader\n")

                for index, header_start in enumerate(header_starts):
                    next_header_start = header_starts[index + 1] if index + 1 < len(header_starts) else len(data)
                    record = index_record_from_mmap(data, header_start, next_header_start)
                    if record.id in seen_ids:
                        raise ValueError(f"duplicate FASTA ID in contig FASTA: {record.id}")
                    seen_ids.add(record.id)
                    out_handle.write(
                        "\t".join(
                            [
                                record.id,
                                str(record.length),
                                str(record.seq_offset if record.seq_offset is not None else -1),
                                str(record.line_bases),
                                str(record.line_width),
                                "1" if record.is_circular else "0",
                                record.header,
                            ]
                        )
                        + "\n"
                    )
                    records_seen += 1
                    total_bases += record.length
                    if progress is not None:
                        progress(0, 1, record.length)

                out_handle.write(f"#bait2contig_record_count={records_seen}\n")
                out_handle.write(f"#bait2contig_total_bases={total_bases}\n")
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise

    os.replace(tmp_path, index_path)
    return records_seen


def find_header_starts(
    data: mmap.mmap,
    *,
    threads: int,
    progress: Optional[FastaProgressCallback] = None,
) -> list[int]:
    """Find FASTA header byte offsets, scanning chunks in parallel."""

    size = len(data)
    if size == 0:
        return []
    ranges = index_scan_ranges(size, threads)
    starts: list[int] = []
    if threads <= 1 or len(ranges) <= 1:
        for start, end in ranges:
            chunk_starts = find_header_starts_in_range(data, start, end)
            starts.extend(chunk_starts)
            if progress is not None:
                progress(end - start, 0, 0)
    else:
        with ThreadPoolExecutor(max_workers=threads, thread_name_prefix="bait2contig-index") as executor:
            futures = {
                executor.submit(find_header_starts_in_range, data, start, end): (start, end)
                for start, end in ranges
            }
            for future in as_completed(futures):
                start, end = futures[future]
                starts.extend(future.result())
                if progress is not None:
                    progress(end - start, 0, 0)
    return sorted(set(starts))


def index_scan_ranges(size: int, threads: int) -> list[tuple[int, int]]:
    chunk_size = 256 * 1024 * 1024
    if size < chunk_size:
        return [(0, size)]
    ranges: list[tuple[int, int]] = []
    for start in range(0, size, chunk_size):
        ranges.append((start, min(size, start + chunk_size)))
    return ranges


def find_header_starts_in_range(data: mmap.mmap, start: int, end: int) -> list[int]:
    starts: list[int] = []
    size = len(data)
    if start == 0 and size and data[0] == ord(">"):
        starts.append(0)
    search = max(0, start - 1)
    stop = min(size, end + 1)
    while True:
        position = data.find(b"\n>", search, stop)
        if position < 0:
            break
        header_start = position + 1
        if start <= header_start < end:
            starts.append(header_start)
        search = header_start + 1
    return starts


def index_record_from_mmap(data: mmap.mmap, header_start: int, record_end: int) -> IndexedFastaRecord:
    header_end = data.find(b"\n", header_start, record_end)
    if header_end < 0:
        header_end = record_end
    header = data[header_start + 1 : header_end].rstrip(b"\r").decode("utf-8").strip()
    if not header:
        raise ValueError("empty FASTA header in contig FASTA")
    record_id = fasta_id_from_header(header)
    seq_offset = min(header_end + 1, record_end)
    first_line_start, line_bases, line_width = first_sequence_line(data, seq_offset, record_end)
    span = max(0, record_end - first_line_start)
    length = sequence_length_from_span(span, line_bases, line_width)
    return IndexedFastaRecord(
        id=record_id,
        header=header,
        length=length,
        is_circular=infer_circular_from_header(header),
        seq_offset=first_line_start,
        line_bases=line_bases,
        line_width=line_width,
    )


def first_sequence_line(data: mmap.mmap, seq_offset: int, record_end: int) -> tuple[int, int, int]:
    line_start = seq_offset
    while line_start < record_end:
        line_end = data.find(b"\n", line_start, record_end)
        has_newline = line_end >= 0
        if not has_newline:
            line_end = record_end
        raw_line = data[line_start:line_end].rstrip(b"\r")
        sequence = raw_line.strip()
        if sequence:
            line_width = line_end - line_start + (1 if has_newline else 0)
            return line_start, len(sequence), line_width
        line_start = line_end + (1 if has_newline else 0)
    return seq_offset, 0, 0


def sequence_length_from_span(span: int, line_bases: int, line_width: int) -> int:
    if span <= 0 or line_bases <= 0:
        return 0
    if line_width <= line_bases:
        return span
    full_lines, remainder = divmod(span, line_width)
    return full_lines * line_bases + min(remainder, line_bases)


class FastaIndex:
    """Read indexed FASTA metadata and fetch plain-FASTA sequences by offset."""

    def __init__(self, source_path: str | Path, index_path: str | Path) -> None:
        self.source_path = Path(source_path)
        self.index_path = Path(index_path)
        self.metadata, self.records = read_index(index_path, load_records=True)
        self.random_access = self.metadata.get("random_access") == "true"

    def close(self) -> None:
        return None

    def __enter__(self) -> "FastaIndex":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def get_info(self, ids: Iterable[str]) -> dict[str, IndexedFastaRecord]:
        wanted = set(ids)
        if not wanted:
            return {}
        return {record_id: record for record_id, record in self.records.items() if record_id in wanted}

    def get_fasta_records(self, ids: Iterable[str]) -> dict[str, FastaRecord]:
        info = self.get_info(ids)
        records: dict[str, FastaRecord] = {}
        if not self.random_access:
            return records
        ordered = sorted(
            (record for record in info.values() if record.seq_offset is not None and record.seq_offset >= 0),
            key=lambda record: int(record.seq_offset or 0),
        )
        with open(self.source_path, "rb") as handle:
            for record in ordered:
                sequence = read_sequence_from_index(handle, record)
                records[record.id] = FastaRecord(
                    id=record.id,
                    header=record.header,
                    sequence=sequence,
                    is_circular=record.is_circular,
                )
        return records


def read_index(index_path: str | Path, *, load_records: bool) -> tuple[dict[str, str], dict[str, IndexedFastaRecord]]:
    metadata: dict[str, str] = {}
    records: dict[str, IndexedFastaRecord] = {}
    with open(index_path, "rt", encoding="utf-8", newline="") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith("#bait2contig_") and "=" in line:
                key, value = line[13:].split("=", 1)
                metadata[key] = value
                continue
            if line.startswith("#"):
                continue
            if not load_records:
                break
            parts = line.split("\t", 6)
            if len(parts) != 7:
                raise ValueError(f"invalid FASTA index record: {line}")
            record_id, length, offset, line_bases, line_width, is_circular, header = parts
            records[record_id] = IndexedFastaRecord(
                id=record_id,
                header=header,
                length=int(length),
                is_circular=is_circular == "1",
                seq_offset=None if int(offset) < 0 else int(offset),
                line_bases=int(line_bases),
                line_width=int(line_width),
            )
    return metadata, records


def read_sequence_from_index(handle, record: IndexedFastaRecord) -> str:
    if record.seq_offset is None:
        return ""
    handle.seek(record.seq_offset)
    remaining = record.length
    parts: list[str] = []
    while remaining > 0:
        raw_line = handle.readline()
        if not raw_line or raw_line.startswith(b">"):
            break
        sequence = raw_line.strip()
        if sequence:
            parts.append(sequence.decode("utf-8"))
            remaining -= len(sequence)
    return "".join(parts)


def read_fasta_subset(path: str | Path, wanted_ids: Iterable[str]) -> dict[str, FastaRecord]:
    """Read selected FASTA records by streaming the input once."""

    wanted = set(wanted_ids)
    if not wanted:
        return {}
    records: dict[str, FastaRecord] = {}
    header: Optional[str] = None
    record_id: Optional[str] = None
    seq_parts: list[str] = []
    keep = False

    def flush() -> None:
        nonlocal header, record_id, seq_parts, keep
        if keep and header is not None and record_id is not None:
            records[record_id] = FastaRecord(
                id=record_id,
                header=header,
                sequence="".join(seq_parts),
                is_circular=infer_circular_from_header(header),
            )
        header = None
        record_id = None
        seq_parts = []
        keep = False

    with open_text(path, "rt") as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                if len(records) == len(wanted):
                    break
                header = line[1:].strip()
                if not header:
                    raise ValueError("empty FASTA header in contig FASTA")
                record_id = fasta_id_from_header(header)
                keep = record_id in wanted
            elif keep:
                seq_parts.append(line.strip())
    flush()
    return records
