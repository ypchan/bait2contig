"""Search workflow for bait2contig."""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

from rich.progress import (
    BarColumn,
    FileSizeColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
    TransferSpeedColumn,
)

from . import __version__
from .io import (
    FastaRecord,
    PafHit,
    copy_or_gzip,
    count_text_lines,
    ensure_parent_dir,
    gzip_output_path,
    parse_paf,
    read_circular_list,
    read_fasta,
    read_lineage,
    write_fasta,
    write_tsv,
)
from .log import (
    DONE_MARKER,
    FAILED_MARKER,
    START_MARKER,
    Logger,
    ResourceMonitor,
    check_resume,
    format_metric,
    now_iso,
)


@dataclass(frozen=True)
class SearchHit:
    """A bait-contig hit after annotation with circularity and lineage."""

    ctg_id: str
    bait_id: str
    identity: float
    aln_length: int
    cov_bait: float
    ctg_len: int
    is_circular: bool
    bait_len: int = 0
    bait_start: int = 0
    bait_end: int = 0
    ctg_start: int = 0
    ctg_end: int = 0
    lineage: Optional[str] = None


class SearchError(RuntimeError):
    """Raised for bait2contig search errors."""


class LoggedSearchError(SearchError):
    """Raised after a search error has already been logged."""


def best_hit_key(hit: SearchHit | PafHit) -> tuple[float, float, int, int]:
    """Return the ranking key used for best-hit selection."""

    return (hit.identity, hit.cov_bait, hit.aln_length, hit.ctg_len)


def filter_hits(
    hits: Iterable[SearchHit],
    *,
    identity: float,
    coverage: float,
    min_aln_length: int,
    terminal_filter: bool = True,
    terminal_tolerance: int = 5,
) -> List[SearchHit]:
    """Filter hits by identity, bait coverage, and alignment length."""

    return [
        hit
        for hit in hits
        if hit.identity >= identity
        and hit.cov_bait >= coverage
        and hit.aln_length >= min_aln_length
        and (not terminal_filter or is_terminal_partial_hit(hit, terminal_tolerance))
    ]


def is_terminal_partial_hit(hit: SearchHit, terminal_tolerance: int = 5) -> bool:
    """Require partial bait alignments to touch sequence ends within tolerance."""

    if hit.bait_len <= 0 or hit.cov_bait >= 1.0:
        return True
    tolerance = max(0, terminal_tolerance)
    bait_at_end = hit.bait_start <= tolerance or (hit.bait_len - hit.bait_end) <= tolerance
    if hit.ctg_len <= 0:
        return bait_at_end
    contig_at_end = hit.ctg_start <= tolerance or (hit.ctg_len - hit.ctg_end) <= tolerance
    return bait_at_end and contig_at_end


def select_best_per_bait(hits: Iterable[SearchHit]) -> List[SearchHit]:
    """Keep one best contig for each bait using the required ranking."""

    best: Dict[str, SearchHit] = {}
    for hit in hits:
        current = best.get(hit.bait_id)
        if current is None or best_hit_key(hit) > best_hit_key(current):
            best[hit.bait_id] = hit
    return list(best.values())


def annotate_paf_hits(
    paf_hits: Iterable[PafHit],
    *,
    contigs: Dict[str, FastaRecord],
    lineage: Optional[Dict[str, str]] = None,
    circular_ids: Optional[set[str]] = None,
) -> List[SearchHit]:
    """Convert parsed PAF records into annotated search hits."""

    annotated: List[SearchHit] = []
    for hit in paf_hits:
        record = contigs.get(hit.ctg_id)
        if circular_ids is not None:
            is_circular = hit.ctg_id in circular_ids
        else:
            is_circular = record.is_circular if record is not None else False
        ctg_len = record.length if record is not None else hit.ctg_len
        annotated.append(
            SearchHit(
                ctg_id=hit.ctg_id,
                bait_id=hit.bait_id,
                identity=hit.identity,
                aln_length=hit.aln_length,
                cov_bait=hit.cov_bait,
                ctg_len=ctg_len,
                is_circular=is_circular,
                bait_len=hit.bait_len,
                bait_start=hit.query_start,
                bait_end=hit.query_end,
                ctg_start=hit.target_start,
                ctg_end=hit.target_end,
                lineage=lineage.get(hit.bait_id, "") if lineage is not None else None,
            )
        )
    return annotated


def search_output_columns(include_lineage: bool) -> List[str]:
    columns = ["ctg_id", "bait_id", "identity", "aln_length", "cov_bait", "ctg_len", "is_circular"]
    if include_lineage:
        columns.append("lineage")
    return columns


def hit_to_row(hit: SearchHit, include_lineage: bool) -> dict[str, object]:
    row: dict[str, object] = {
        "ctg_id": hit.ctg_id,
        "bait_id": hit.bait_id,
        "identity": f"{hit.identity:.6f}",
        "aln_length": hit.aln_length,
        "cov_bait": f"{hit.cov_bait:.6f}",
        "ctg_len": hit.ctg_len,
        "is_circular": str(hit.is_circular),
    }
    if include_lineage:
        row["lineage"] = hit.lineage or ""
    return row


def write_hits_tsv(path: str | Path, hits: Iterable[SearchHit], *, include_lineage: bool) -> None:
    columns = search_output_columns(include_lineage)
    write_tsv(path, columns, [hit_to_row(hit, include_lineage) for hit in hits])


def select_extraction_hits(
    hits: Sequence[SearchHit],
    *,
    mode: str,
    identity: float,
    coverage: float,
    min_aln_length: int,
    terminal_filter: bool = True,
    terminal_tolerance: int = 5,
) -> List[SearchHit]:
    selected = filter_hits(
        hits,
        identity=identity,
        coverage=coverage,
        min_aln_length=min_aln_length,
        terminal_filter=terminal_filter,
        terminal_tolerance=terminal_tolerance,
    )
    if mode == "best":
        selected = select_best_per_bait(selected)
    elif mode == "circular":
        selected = [hit for hit in selected if hit.is_circular]
    elif mode == "non-circular":
        selected = [hit for hit in selected if not hit.is_circular]
    elif mode != "all":
        raise SearchError(f"unsupported extract mode: {mode}")
    return selected


def renamed_header(hit: SearchHit, *, include_lineage: bool = False) -> str:
    header = (
        f"{hit.ctg_id}|bait={hit.bait_id}|identity={hit.identity:.6f}|"
        f"cov_bait={hit.cov_bait:.6f}|aln_length={hit.aln_length}|"
        f"ctg_len={hit.ctg_len}|circular={hit.is_circular}"
    )
    if include_lineage:
        header += f"|lineage={hit.lineage or ''}"
    return header


def extract_contigs(
    hits: Sequence[SearchHit],
    contigs: Dict[str, FastaRecord],
    output_path: str | Path,
    *,
    rename: bool,
    include_lineage: bool,
    dedup: bool,
    logger: Optional[Logger] = None,
) -> int:
    """Write matched contigs to FASTA and return the number written."""

    if include_lineage and not rename:
        raise SearchError("--extract-include-lineage requires --extract-rename.")
    if not dedup and not rename and logger is not None:
        logger.warn("Duplicate FASTA IDs may be written because --no-extract-dedup was used without --extract-rename.")

    seen: set[str] = set()
    fasta_records: List[tuple[str, str]] = []
    for hit in hits:
        record = contigs.get(hit.ctg_id)
        if record is None:
            if logger is not None:
                logger.warn(f"matched contig is absent from contig FASTA and cannot be extracted: {hit.ctg_id}")
            continue
        if dedup and hit.ctg_id in seen:
            continue
        seen.add(hit.ctg_id)
        header = renamed_header(hit, include_lineage=include_lineage) if rename else record.header
        fasta_records.append((header, record.sequence))
    write_fasta(fasta_records, output_path)
    return len(fasta_records)


def final_paf_path(actual_out: str | Path, gzip_enabled: bool) -> str:
    """Return the final kept PAF path derived from the actual output path."""

    path = str(actual_out)
    if path.endswith(".gz"):
        path = path[:-3]
    if path.endswith(".tsv"):
        path = path[:-4]
    path = f"{path}.paf"
    if gzip_enabled:
        path = gzip_output_path(path, True)
    return path


def resolve_minimap2(path: str) -> Optional[str]:
    if os.sep in path or (os.altsep and os.altsep in path):
        return path if Path(path).exists() else None
    return shutil.which(path)


def minimap2_version(executable: Optional[str]) -> str:
    if executable is None:
        return "NA"
    try:
        completed = subprocess.run(
            [executable, "--version"],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except Exception:
        return "NA"
    version = (completed.stdout or completed.stderr).strip().splitlines()
    return version[0] if version else "NA"


def abs_or_na(path: Optional[str]) -> str:
    return str(Path(path).resolve()) if path else "NA"


def optional_value(value: object) -> object:
    return value if value is not None else "NA"


def require_existing_file(path: str, option: str) -> None:
    file_path = Path(path)
    if not file_path.exists():
        raise SearchError(f"{option} file not found: {path}")
    if not file_path.is_file():
        raise SearchError(f"{option} must be a file: {path}")


def validate_output_path(path: Optional[str], option: str) -> None:
    if path is None:
        return
    output_path = Path(path)
    if output_path.exists() and output_path.is_dir():
        raise SearchError(f"{option} must be a file path, not a directory: {path}")
    parent = output_path.parent
    if str(parent) and parent.exists() and not parent.is_dir():
        raise SearchError(f"{option} parent path is not a directory: {parent}")


def build_search_params(args, actual_out: str, actual_extract: Optional[str], minimap2_version_value: str) -> dict[str, object]:
    return {
        "command": "search",
        "status": "running",
        "started_at": now_iso(),
        "bait2contig_version": __version__,
        "python_version": ".".join(map(str, tuple(__import__("sys").version_info[:3]))),
        "out": str(Path(actual_out).resolve()),
        "gzip": bool(args.gzip),
        "contigs": str(Path(args.contigs).resolve()),
        "bait": str(Path(args.bait).resolve()),
        "lineage": abs_or_na(args.lineage),
        "circular_list": abs_or_na(args.circular_list),
        "identity": args.identity,
        "coverage": args.coverage,
        "min_aln_length": args.min_aln_length,
        "best_only": bool(args.best_only),
        "terminal_filter": bool(args.terminal_filter),
        "terminal_tolerance": args.terminal_tolerance,
        "preset": args.preset,
        "threads": args.threads,
        "minimap2_version": minimap2_version_value,
        "extract_contigs": abs_or_na(actual_extract),
        "extract_mode": args.extract_mode,
        "extract_min_identity": args.extract_min_identity,
        "extract_min_coverage": args.extract_min_coverage,
        "extract_min_aln_length": args.extract_min_aln_length,
        "extract_rename": bool(args.extract_rename),
        "extract_dedup": bool(args.extract_dedup),
        "extract_include_lineage": bool(args.extract_include_lineage),
    }


def search_resume_params(start_params: dict[str, object]) -> dict[str, object]:
    keys = [
        "out",
        "gzip",
        "contigs",
        "bait",
        "lineage",
        "circular_list",
        "identity",
        "coverage",
        "min_aln_length",
        "best_only",
        "terminal_filter",
        "terminal_tolerance",
        "preset",
        "minimap2_version",
        "extract_contigs",
        "extract_mode",
        "extract_min_identity",
        "extract_min_coverage",
        "extract_min_aln_length",
        "extract_rename",
        "extract_dedup",
        "extract_include_lineage",
    ]
    return {key: start_params[key] for key in keys}


def validate_search_args(args) -> None:
    if args.resume and args.rerun:
        raise SearchError("--resume and --rerun cannot be used together.")
    if args.extract_include_lineage and not args.extract_rename:
        raise SearchError("--extract-include-lineage requires --extract-rename.")
    require_existing_file(args.contigs, "--contigs")
    require_existing_file(args.bait, "--bait")
    if args.lineage:
        require_existing_file(args.lineage, "--lineage")
    if args.circular_list:
        require_existing_file(args.circular_list, "--circular-list")
    validate_output_path(args.out, "--out")
    validate_output_path(args.extract_contigs, "--extract-contigs")
    validate_output_path(args.log, "--log")
    if args.tmp_dir and Path(args.tmp_dir).exists() and not Path(args.tmp_dir).is_dir():
        raise SearchError(f"--tmp-dir must be a directory: {args.tmp_dir}")
    for label in ("identity", "coverage"):
        value = getattr(args, label)
        if value < 0 or value > 1:
            raise SearchError(f"--{label.replace('_', '-')} must be between 0 and 1")
    for label in ("extract_min_identity", "extract_min_coverage"):
        value = getattr(args, label)
        if value is not None and (value < 0 or value > 1):
            raise SearchError(f"--{label.replace('_', '-')} must be between 0 and 1")
    if args.min_aln_length < 0:
        raise SearchError("--min-aln-length must be at least 0")
    if args.terminal_tolerance < 0:
        raise SearchError("--terminal-tolerance must be at least 0")
    if args.extract_min_aln_length is not None and args.extract_min_aln_length < 0:
        raise SearchError("--extract-min-aln-length must be at least 0")
    if args.threads < 1:
        raise SearchError("--threads must be at least 1")
    if args.monitor_interval < 1:
        raise SearchError("--monitor-interval must be at least 1")


def check_existing_outputs(paths: Iterable[Optional[str]], *, resume: bool, rerun: bool, force: bool) -> None:
    if resume or rerun or force:
        return
    for path in paths:
        if path and Path(path).exists():
            raise SearchError(
                f"output already exists: {path}\n"
                "Use --resume to reuse a completed run, --rerun to recompute, or --force to overwrite."
            )


def run_minimap2(
    *,
    executable: str,
    preset: str,
    threads: int,
    contigs: str,
    bait: str,
    tmp_paf: str,
    logger: Logger,
    monitor: ResourceMonitor,
) -> None:
    command = [executable, "-x", preset, "-t", str(threads), contigs, bait]
    logger.info("running minimap2")
    with open(tmp_paf, "wt", encoding="utf-8", newline="") as paf_handle:
        process = subprocess.Popen(
            command,
            stdout=paf_handle,
            stderr=subprocess.PIPE,
            text=True,
        )
        monitor.set_child_pid(process.pid)
        _, stderr = process.communicate()
    if stderr:
        for line in stderr.splitlines():
            if line.strip() and logger.verbose:
                logger.info(f"minimap2: {line.strip()}")
    if process.returncode != 0:
        raise SearchError(f"minimap2 failed with exit code {process.returncode}")


def make_tmp_paf(tmp_dir: str | Path) -> str:
    Path(tmp_dir).mkdir(parents=True, exist_ok=True)
    fd, path = tempfile.mkstemp(prefix="bait2contig.", suffix=".paf", dir=str(tmp_dir))
    os.close(fd)
    return path


def fasta_progress_total(path: str | Path) -> Optional[int]:
    text_path = str(path)
    if text_path.endswith(".gz"):
        return None
    try:
        return Path(path).stat().st_size
    except OSError:
        return None


def format_bases(count: int) -> str:
    if count >= 1_000_000_000:
        return f"{count / 1_000_000_000:.2f} Gb"
    if count >= 1_000_000:
        return f"{count / 1_000_000:.2f} Mb"
    if count >= 1_000:
        return f"{count / 1_000:.2f} kb"
    return f"{count} bp"


class FastaProgressTracker:
    """Throttle FASTA parser progress updates to keep large inputs responsive."""

    def __init__(self, progress: Progress, task_id: int, total: Optional[int]) -> None:
        self.progress = progress
        self.task_id = task_id
        self.total = total
        self.bytes_seen = 0
        self.records_seen = 0
        self.bases_seen = 0
        self._last_update = time.time()
        self._last_update_bytes = 0

    def __call__(self, bytes_delta: int, records_delta: int, bases_delta: int) -> None:
        self.bytes_seen += bytes_delta
        self.records_seen += records_delta
        self.bases_seen += bases_delta
        now = time.time()
        if self.bytes_seen - self._last_update_bytes >= 1_048_576 or now - self._last_update >= 0.25:
            self.flush()

    def flush(self) -> None:
        completed = self.bytes_seen
        if self.total is not None:
            completed = min(completed, self.total)
        self.progress.update(
            self.task_id,
            completed=completed,
            records=f"{self.records_seen:,}",
            bases=format_bases(self.bases_seen),
        )
        self._last_update = time.time()
        self._last_update_bytes = self.bytes_seen


def make_fasta_progress(logger: Logger) -> Progress:
    return Progress(
        SpinnerColumn(),
        TextColumn("[bold]{task.description}"),
        BarColumn(),
        FileSizeColumn(),
        TransferSpeedColumn(),
        TimeElapsedColumn(),
        TextColumn("{task.fields[records]} seq"),
        TextColumn("{task.fields[bases]}"),
        console=logger.console,
        transient=True,
        disable=not logger.progress_enabled,
    )


def load_fasta_inputs(
    args,
    *,
    logger: Logger,
    monitor: ResourceMonitor,
) -> tuple[Dict[str, FastaRecord], Dict[str, FastaRecord]]:
    """Read bait and contig FASTA inputs with progress and parallel file loading."""

    monitor.set_stage("loading_fasta_inputs")
    logger.info("loading FASTA inputs")
    inputs = {
        "bait": args.bait,
        "contigs": args.contigs,
    }
    results: dict[str, Dict[str, FastaRecord]] = {}

    with make_fasta_progress(logger) as progress:
        trackers: dict[str, FastaProgressTracker] = {}
        for name, path in inputs.items():
            total = fasta_progress_total(path)
            task_id = progress.add_task(
                f"{name} FASTA",
                total=total,
                records="0",
                bases="0 bp",
            )
            trackers[name] = FastaProgressTracker(progress, task_id, total)

        def read_one(name: str, path: str) -> tuple[str, Dict[str, FastaRecord]]:
            tracker = trackers[name]
            try:
                return name, read_fasta(path, progress=tracker)
            finally:
                tracker.flush()

        with ThreadPoolExecutor(max_workers=2, thread_name_prefix="bait2contig-fasta") as executor:
            futures = [executor.submit(read_one, name, path) for name, path in inputs.items()]
            for future in as_completed(futures):
                name, records = future.result()
                results[name] = records

    return results["bait"], results["contigs"]


def run_search(args) -> int:
    """Run bait2contig search from parsed argparse arguments."""

    validate_search_args(args)
    args.extract_min_identity = args.identity if args.extract_min_identity is None else args.extract_min_identity
    args.extract_min_coverage = args.coverage if args.extract_min_coverage is None else args.extract_min_coverage
    args.extract_min_aln_length = (
        args.min_aln_length if args.extract_min_aln_length is None else args.extract_min_aln_length
    )

    actual_out = gzip_output_path(args.out, args.gzip)
    actual_extract = gzip_output_path(args.extract_contigs, args.gzip) if args.extract_contigs else None
    log_path = args.log or f"{actual_out}.log"
    tmp_dir = args.tmp_dir or str(Path(actual_out).resolve().parent)

    executable_for_resume = resolve_minimap2(args.minimap2)
    minimap2_version_value = minimap2_version(executable_for_resume)
    initial_params = build_search_params(args, actual_out, actual_extract, minimap2_version_value)

    if args.resume:
        with Logger(log_path, no_color=args.no_color, quiet=args.quiet, verbose=args.verbose) as logger:
            result = check_resume(
                log_path=log_path,
                command="search",
                output_path=actual_out,
                expected_params=search_resume_params(initial_params),
            )
            if result.ok:
                logger.info("Resume enabled: previous successful run found in log.")
                logger.info(f"Output verified: {Path(actual_out).name}")
                logger.info("Skipping search.")
                return 0
            logger.warn(f"Resume log found but {result.reason}")
            logger.info("Re-running search.")
    else:
        check_existing_outputs(
            [actual_out, actual_extract],
            resume=args.resume,
            rerun=args.rerun,
            force=args.force,
        )

    ensure_parent_dir(actual_out)
    if actual_extract:
        ensure_parent_dir(actual_extract)
    logger = Logger(log_path, no_color=args.no_color, quiet=args.quiet, verbose=args.verbose)
    monitor = ResourceMonitor(args.monitor_interval, logger)
    monitor.start()
    start_time = time.time()
    raw_alignment_count = 0
    kept_alignment_count = 0
    extracted_count = 0
    extract_output_size = 0
    tmp_paf: Optional[str] = None
    start_params = initial_params
    try:
        logger.info(f"bait2contig version: {__version__}")
        logger.info("command: search")
        logger.info(f"output: {actual_out}")
        logger.info(f"log: {log_path}")
        logger.info(f"resume: {'enabled' if args.resume else 'disabled'}")
        logger.info(f"rerun: {'enabled' if args.rerun else 'disabled'}")
        logger.info(f"resource monitor backend: {monitor.backend}")
        logger.marker(START_MARKER, start_params)

        executable = executable_for_resume
        if executable is None:
            raise SearchError("minimap2 was not found. Please install minimap2 or provide its path with --minimap2.")

        bait, contigs = load_fasta_inputs(args, logger=logger, monitor=monitor)
        logger.info(f"loaded bait sequences: {len(bait)}")
        logger.info(f"loaded contigs: {len(contigs)}")

        monitor.set_stage("loading_annotations")
        lineage = read_lineage(args.lineage) if args.lineage else None
        if lineage is not None:
            extra_lineage_ids = sorted(set(lineage) - set(bait))
            for bait_id in extra_lineage_ids:
                logger.warn(f"lineage contains bait_id not present in bait FASTA: {bait_id}")
        circular_ids = read_circular_list(args.circular_list) if args.circular_list else None

        tmp_paf = make_tmp_paf(tmp_dir)
        monitor.set_stage("running_minimap2")
        run_minimap2(
            executable=executable,
            preset=args.preset,
            threads=args.threads,
            contigs=args.contigs,
            bait=args.bait,
            tmp_paf=tmp_paf,
            logger=logger,
            monitor=monitor,
        )

        monitor.set_stage("parsing_paf")
        logger.info("parsing PAF")
        paf_hits = parse_paf(tmp_paf)
        raw_alignment_count = len(paf_hits)
        monitor.set_stage("filtering_hits")
        all_hits = annotate_paf_hits(paf_hits, contigs=contigs, lineage=lineage, circular_ids=circular_ids)
        kept_hits = filter_hits(
            all_hits,
            identity=args.identity,
            coverage=args.coverage,
            min_aln_length=args.min_aln_length,
            terminal_filter=args.terminal_filter,
            terminal_tolerance=args.terminal_tolerance,
        )
        if args.best_only:
            kept_hits = select_best_per_bait(kept_hits)
        kept_alignment_count = len(kept_hits)
        logger.info(f"raw alignments: {raw_alignment_count}")
        logger.info(f"kept alignments: {kept_alignment_count}")

        monitor.set_stage("writing_hits")
        write_hits_tsv(actual_out, kept_hits, include_lineage=lineage is not None)

        if args.keep_paf:
            paf_output = final_paf_path(actual_out, args.gzip)
            copy_or_gzip(tmp_paf, paf_output, args.gzip)
            tmp_paf = None
            logger.info(f"kept PAF: {paf_output}")
        elif tmp_paf:
            Path(tmp_paf).unlink(missing_ok=True)
            tmp_paf = None

        if actual_extract:
            monitor.set_stage("extracting_contigs")
            extraction_hits = select_extraction_hits(
                all_hits,
                mode=args.extract_mode,
                identity=args.extract_min_identity,
                coverage=args.extract_min_coverage,
                min_aln_length=args.extract_min_aln_length,
                terminal_filter=args.terminal_filter,
                terminal_tolerance=args.terminal_tolerance,
            )
            extracted_count = extract_contigs(
                extraction_hits,
                contigs,
                actual_extract,
                rename=args.extract_rename,
                include_lineage=args.extract_include_lineage,
                dedup=args.extract_dedup,
                logger=logger,
            )
            extract_output_size = Path(actual_extract).stat().st_size if Path(actual_extract).exists() else 0

        stats = monitor.stop()
        runtime = time.time() - start_time
        output_size = Path(actual_out).stat().st_size if Path(actual_out).exists() else 0
        output_lines = count_text_lines(actual_out) if Path(actual_out).exists() else 0
        done_values = {
            "command": "search",
            "status": "success",
            "finished_at": now_iso(),
            "runtime_seconds": f"{runtime:.2f}",
            "exit_code": 0,
            "output": str(Path(actual_out).resolve()),
            "output_size": output_size,
            "output_lines": output_lines,
            "raw_alignment_count": raw_alignment_count,
            "kept_alignment_count": kept_alignment_count,
            "unique_bait_hit_count": len({hit.bait_id for hit in kept_hits}),
            "unique_contig_hit_count": len({hit.ctg_id for hit in kept_hits}),
            "extracted_contig_count": extracted_count,
            "extract_output": str(Path(actual_extract).resolve()) if actual_extract else "NA",
            "extract_output_size": extract_output_size,
            "peak_rss_mb": format_metric(stats["peak_rss_mb"], digits=2),
            "mean_cpu_percent": format_metric(stats["mean_cpu_percent"], digits=1),
            "max_cpu_percent": format_metric(stats["max_cpu_percent"], digits=1),
        }
        logger.marker(DONE_MARKER, done_values)
        logger.done(f"wrote output: {actual_out}")
        logger.done("bait2contig search completed successfully")
        logger.done(f"runtime: {runtime:.2f} sec")
        logger.done(f"peak RSS: {format_metric(stats['peak_rss_mb'], digits=2)} MB")
        logger.done(f"mean CPU: {format_metric(stats['mean_cpu_percent'], digits=1)}%")
        logger.close()
        return 0
    except Exception as exc:
        if tmp_paf:
            Path(tmp_paf).unlink(missing_ok=True)
        stats = monitor.stop()
        runtime = time.time() - start_time
        message = str(exc)
        logger.error(f"ERROR: {message}")
        logger.marker(
            FAILED_MARKER,
            {
                "command": "search",
                "status": "failed",
                "finished_at": now_iso(),
                "runtime_seconds": f"{runtime:.2f}",
                "exit_code": 1,
                "error_type": type(exc).__name__,
                "error": message,
                "peak_rss_mb": format_metric(stats["peak_rss_mb"], digits=2),
                "mean_cpu_percent": format_metric(stats["mean_cpu_percent"], digits=1),
                "max_cpu_percent": format_metric(stats["max_cpu_percent"], digits=1),
            },
        )
        logger.close()
        raise LoggedSearchError(message) from exc
