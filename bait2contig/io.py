"""Input and output helpers for bait2contig."""

from __future__ import annotations

import csv
import gzip
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, Iterator, List, Optional, Sequence, TextIO, Tuple


@dataclass(frozen=True)
class FastaRecord:
    """A parsed FASTA record."""

    id: str
    header: str
    sequence: str
    is_circular: bool = False

    @property
    def length(self) -> int:
        return len(self.sequence)


@dataclass(frozen=True)
class PafHit:
    """A parsed PAF alignment normalized to bait2contig semantic fields."""

    bait_id: str
    bait_len: int
    bait_start: int
    bait_end: int
    strand: str
    ctg_id: str
    ctg_len: int
    ctg_start: int
    ctg_end: int
    residue_matches: int
    aln_length: int
    mapping_quality: int
    identity: float
    cov_bait: float


def open_text(path: str | Path, mode: str = "rt") -> TextIO:
    """Open plain-text or gzip-compressed text files with UTF-8 encoding."""

    path = str(path)
    if "b" in mode:
        raise ValueError("open_text only supports text modes")
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")


def gzip_output_path(path: str | Path, gzip_enabled: bool) -> str:
    """Return the actual output path after applying a --gzip request."""

    text_path = str(path)
    if gzip_enabled and not text_path.endswith(".gz"):
        return f"{text_path}.gz"
    return text_path


def ensure_parent_dir(path: str | Path) -> None:
    parent = Path(path).parent
    if str(parent) and not parent.exists():
        parent.mkdir(parents=True, exist_ok=True)


def fasta_id_from_header(header: str) -> str:
    """Return the first whitespace-delimited token from a FASTA header."""

    stripped = header.strip()
    if not stripped:
        raise ValueError("FASTA header is empty")
    return stripped.split()[0]


def infer_circular_from_header(header: str) -> bool:
    """Infer circularity from supported FASTA header keywords."""

    normalized = header.lower().replace(";", " ").replace(",", " ")
    tokens = [token.strip() for token in normalized.split() if token.strip()]
    circular_tokens = {"circular", "circular=true", "is_circular=true", "circ=true"}
    return any(token in circular_tokens for token in tokens)


FastaProgressCallback = Callable[[int, int, int], None]


def read_fasta(path: str | Path, progress: Optional[FastaProgressCallback] = None) -> Dict[str, FastaRecord]:
    """Read a plain or gzip FASTA file into records keyed by first-token ID."""

    records: Dict[str, FastaRecord] = {}
    header: Optional[str] = None
    seq_parts: List[str] = []

    def flush() -> None:
        nonlocal header, seq_parts
        if header is None:
            return
        record_id = fasta_id_from_header(header)
        if record_id in records:
            raise ValueError(f"duplicate FASTA ID: {record_id}")
        records[record_id] = FastaRecord(
            id=record_id,
            header=header,
            sequence="".join(seq_parts),
            is_circular=infer_circular_from_header(header),
        )
        if progress is not None:
            progress(0, 1, 0)
        header = None
        seq_parts = []

    with open_text(path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if progress is not None:
                progress(len(line), 0, 0)
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"empty FASTA header at line {line_number}")
            else:
                if header is None:
                    raise ValueError(f"sequence found before FASTA header at line {line_number}")
                sequence = line.strip()
                seq_parts.append(sequence)
                if progress is not None:
                    progress(0, 0, len(sequence))
    flush()
    return records


def read_fasta_ids(path: str | Path, progress: Optional[FastaProgressCallback] = None) -> set[str]:
    """Read only FASTA record IDs without storing sequences."""

    ids: set[str] = set()
    with open_text(path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if progress is not None:
                progress(len(line), 0, 0)
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            if not header:
                raise ValueError(f"empty FASTA header at line {line_number}")
            record_id = fasta_id_from_header(header)
            if record_id in ids:
                raise ValueError(f"duplicate FASTA ID: {record_id}")
            ids.add(record_id)
            if progress is not None:
                progress(0, 1, 0)
    return ids


def write_fasta(records: Iterable[Tuple[str, str]], path: str | Path, line_width: int = 80) -> None:
    """Write FASTA records as (header, sequence) tuples."""

    ensure_parent_dir(path)
    with open_text(path, "wt") as handle:
        for header, sequence in records:
            clean_header = header[1:] if header.startswith(">") else header
            handle.write(f">{clean_header}\n")
            if not sequence:
                handle.write("\n")
                continue
            for start in range(0, len(sequence), line_width):
                handle.write(f"{sequence[start:start + line_width]}\n")


def read_lineage(path: str | Path) -> Dict[str, str]:
    """Read a two-column bait lineage TSV with or without a header."""

    lineage: Dict[str, str] = {}
    first_data_row = True
    with open_text(path, "rt") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n\r")
            if not line:
                continue
            parts = line.split("\t", 1)
            if len(parts) < 2:
                raise ValueError(f"lineage TSV line {line_number} must contain bait_id and lineage")
            bait_id, value = parts[0].strip(), parts[1].strip()
            if first_data_row and bait_id.lower() == "bait_id" and value.lower() == "lineage":
                first_data_row = False
                continue
            first_data_row = False
            if not bait_id:
                raise ValueError(f"lineage TSV line {line_number} has an empty bait_id")
            if bait_id in lineage:
                raise ValueError(f"duplicate bait_id in lineage file: {bait_id}")
            lineage[bait_id] = value
    return lineage


def read_circular_list(path: str | Path) -> set[str]:
    """Read a plain or gzip list of circular contig IDs."""

    contigs: set[str] = set()
    with open_text(path, "rt") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line:
                contigs.add(line.split()[0])
    return contigs


def parse_paf_line(line: str, line_number: int = 0, *, query_is_bait: bool = False) -> PafHit:
    """Parse one minimap2 PAF line and compute identity and bait coverage.

    By default, bait2contig expects minimap2 to run with contigs as the query
    and bait/reference sequences as the target. Set query_is_bait=True only for
    legacy PAF produced in the older bait-query orientation.
    """

    fields = line.rstrip("\n\r").split("\t")
    if len(fields) < 12:
        where = f" at line {line_number}" if line_number else ""
        raise ValueError(f"PAF record{where} has fewer than 12 columns")
    query_id = fields[0]
    query_len = int(fields[1])
    query_start = int(fields[2])
    query_end = int(fields[3])
    strand = fields[4]
    target_id = fields[5]
    target_len = int(fields[6])
    target_start = int(fields[7])
    target_end = int(fields[8])
    residue_matches = int(fields[9])
    aln_length = int(fields[10])
    mapping_quality = int(fields[11])
    if query_is_bait:
        bait_id = query_id
        bait_len = query_len
        bait_start = query_start
        bait_end = query_end
        ctg_id = target_id
        ctg_len = target_len
        ctg_start = target_start
        ctg_end = target_end
    else:
        ctg_id = query_id
        ctg_len = query_len
        ctg_start = query_start
        ctg_end = query_end
        bait_id = target_id
        bait_len = target_len
        bait_start = target_start
        bait_end = target_end
    identity = residue_matches / aln_length if aln_length else 0.0
    cov_bait = (bait_end - bait_start) / bait_len if bait_len else 0.0
    return PafHit(
        bait_id=bait_id,
        bait_len=bait_len,
        bait_start=bait_start,
        bait_end=bait_end,
        strand=strand,
        ctg_id=ctg_id,
        ctg_len=ctg_len,
        ctg_start=ctg_start,
        ctg_end=ctg_end,
        residue_matches=residue_matches,
        aln_length=aln_length,
        mapping_quality=mapping_quality,
        identity=identity,
        cov_bait=cov_bait,
    )


def parse_paf(path: str | Path, *, query_is_bait: bool = False) -> List[PafHit]:
    """Parse a plain-text PAF file."""

    hits: List[PafHit] = []
    with open_text(path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            hits.append(parse_paf_line(line, line_number, query_is_bait=query_is_bait))
    return hits


def read_tsv(path: str | Path) -> Tuple[List[str], List[dict[str, str]]]:
    """Read a tab-delimited file into a header and string row dictionaries."""

    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return [], []
        rows = [dict(row) for row in reader]
        return list(reader.fieldnames), rows


def write_tsv(path: str | Path, columns: Sequence[str], rows: Iterable[dict[str, object]]) -> None:
    """Write rows to a tab-delimited file."""

    ensure_parent_dir(path)
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def count_text_lines(path: str | Path) -> int:
    """Count lines in a plain or gzip text file."""

    with open_text(path, "rt") as handle:
        return sum(1 for _ in handle)


def copy_or_gzip(src: str | Path, dst: str | Path, gzip_enabled: bool) -> None:
    """Move or gzip a temporary plain-text file to its final destination."""

    ensure_parent_dir(dst)
    if gzip_enabled:
        with open(src, "rb") as in_handle, gzip.open(dst, "wb") as out_handle:
            shutil.copyfileobj(in_handle, out_handle)
        Path(src).unlink(missing_ok=True)
    else:
        shutil.move(str(src), str(dst))
