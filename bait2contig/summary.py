"""Summarize bait2contig hit-level TSV files."""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

from . import __version__
from .core import best_hit_key
from .io import count_text_lines, ensure_parent_dir, gzip_output_path, read_tsv, write_tsv
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
class SummaryHit:
    ctg_id: str
    bait_id: str
    identity: float
    aln_length: int
    cov_bait: float
    ctg_len: int
    is_circular: bool
    lineage: Optional[str] = None


class SummaryError(RuntimeError):
    """Raised for bait2contig summarize errors."""


class LoggedSummaryError(SummaryError):
    """Raised after a summarize error has already been logged."""


def parse_bool(value: str) -> bool:
    return value.strip().lower() == "true"


def summary_best_key(hit: SummaryHit) -> tuple[float, float, int, int]:
    return (hit.identity, hit.cov_bait, hit.aln_length, hit.ctg_len)


def read_hit_rows(path: str | Path) -> tuple[List[str], List[SummaryHit]]:
    columns, rows = read_tsv(path)
    if not columns:
        return [], []
    required = {"ctg_id", "bait_id", "identity", "aln_length", "cov_bait", "ctg_len", "is_circular"}
    missing = required - set(columns)
    if missing:
        raise SummaryError(f"hits TSV is missing required columns: {', '.join(sorted(missing))}")
    hits: List[SummaryHit] = []
    for row in rows:
        hits.append(
            SummaryHit(
                ctg_id=row["ctg_id"],
                bait_id=row["bait_id"],
                identity=float(row["identity"]),
                aln_length=int(row["aln_length"]),
                cov_bait=float(row["cov_bait"]),
                ctg_len=int(row["ctg_len"]),
                is_circular=parse_bool(row["is_circular"]),
                lineage=row.get("lineage") if "lineage" in columns else None,
            )
        )
    return columns, hits


def filter_summary_hits(
    hits: Iterable[SummaryHit],
    *,
    min_identity: Optional[float] = None,
    min_coverage: Optional[float] = None,
    min_aln_length: Optional[int] = None,
) -> List[SummaryHit]:
    kept: List[SummaryHit] = []
    for hit in hits:
        if min_identity is not None and hit.identity < min_identity:
            continue
        if min_coverage is not None and hit.cov_bait < min_coverage:
            continue
        if min_aln_length is not None and hit.aln_length < min_aln_length:
            continue
        kept.append(hit)
    return kept


def deduplicate_bait_contig_hits(hits: Iterable[SummaryHit]) -> List[SummaryHit]:
    best: Dict[tuple[str, str], SummaryHit] = {}
    order: List[tuple[str, str]] = []
    for hit in hits:
        key = (hit.bait_id, hit.ctg_id)
        if key not in best:
            order.append(key)
            best[key] = hit
        elif summary_best_key(hit) > summary_best_key(best[key]):
            best[key] = hit
    return [best[key] for key in order]


def summary_columns(*, include_lineage: bool, best_hit: bool, include_contigs: bool) -> List[str]:
    columns = ["bait_id"]
    if include_lineage:
        columns.append("lineage")
    columns.extend(
        [
            "contig_count",
            "circular_contig_count",
            "total_ctg_len",
            "mean_ctg_len",
            "max_ctg_len",
            "mean_identity",
            "max_identity",
            "mean_cov_bait",
            "max_cov_bait",
            "mean_aln_length",
            "max_aln_length",
        ]
    )
    if best_hit:
        columns.extend(
            [
                "best_ctg_id",
                "best_identity",
                "best_cov_bait",
                "best_aln_length",
                "best_ctg_len",
                "best_is_circular",
            ]
        )
    if include_contigs:
        columns.append("contigs")
    return columns


def mean(values: Sequence[float | int]) -> float:
    return float(sum(values)) / len(values)


def summarize_hits(
    hits: Iterable[SummaryHit],
    *,
    include_lineage: bool,
    best_hit: bool = False,
    include_contigs: bool = False,
    contig_sep: str = ",",
) -> List[dict[str, object]]:
    deduped = deduplicate_bait_contig_hits(hits)
    grouped: Dict[str, List[SummaryHit]] = {}
    bait_order: List[str] = []
    for hit in deduped:
        if hit.bait_id not in grouped:
            grouped[hit.bait_id] = []
            bait_order.append(hit.bait_id)
        grouped[hit.bait_id].append(hit)

    rows: List[dict[str, object]] = []
    for bait_id in bait_order:
        bait_hits = sorted(grouped[bait_id], key=summary_best_key, reverse=True)
        ctg_lengths = [hit.ctg_len for hit in bait_hits]
        identities = [hit.identity for hit in bait_hits]
        coverages = [hit.cov_bait for hit in bait_hits]
        aln_lengths = [hit.aln_length for hit in bait_hits]
        row: dict[str, object] = {
            "bait_id": bait_id,
            "contig_count": len(bait_hits),
            "circular_contig_count": sum(1 for hit in bait_hits if hit.is_circular),
            "total_ctg_len": sum(ctg_lengths),
            "mean_ctg_len": f"{mean(ctg_lengths):.2f}",
            "max_ctg_len": max(ctg_lengths),
            "mean_identity": f"{mean(identities):.6f}",
            "max_identity": f"{max(identities):.6f}",
            "mean_cov_bait": f"{mean(coverages):.6f}",
            "max_cov_bait": f"{max(coverages):.6f}",
            "mean_aln_length": f"{mean(aln_lengths):.2f}",
            "max_aln_length": max(aln_lengths),
        }
        if include_lineage:
            row["lineage"] = next((hit.lineage or "" for hit in bait_hits if hit.lineage is not None), "")
        if best_hit:
            best = bait_hits[0]
            row.update(
                {
                    "best_ctg_id": best.ctg_id,
                    "best_identity": f"{best.identity:.6f}",
                    "best_cov_bait": f"{best.cov_bait:.6f}",
                    "best_aln_length": best.aln_length,
                    "best_ctg_len": best.ctg_len,
                    "best_is_circular": str(best.is_circular),
                }
            )
        if include_contigs:
            row["contigs"] = contig_sep.join(hit.ctg_id for hit in bait_hits)
        rows.append(row)
    return rows


def summarize_file(
    *,
    hits_path: str | Path,
    out_path: str | Path,
    min_identity: Optional[float] = None,
    min_coverage: Optional[float] = None,
    min_aln_length: Optional[int] = None,
    best_hit: bool = False,
    include_contigs: bool = False,
    contig_sep: str = ",",
) -> dict[str, int]:
    input_columns, hits = read_hit_rows(hits_path)
    include_lineage = "lineage" in input_columns
    filtered = filter_summary_hits(
        hits,
        min_identity=min_identity,
        min_coverage=min_coverage,
        min_aln_length=min_aln_length,
    )
    deduped = deduplicate_bait_contig_hits(filtered)
    rows = summarize_hits(
        deduped,
        include_lineage=include_lineage,
        best_hit=best_hit,
        include_contigs=include_contigs,
        contig_sep=contig_sep,
    )
    columns = summary_columns(include_lineage=include_lineage, best_hit=best_hit, include_contigs=include_contigs)
    write_tsv(out_path, columns, rows)
    return {
        "input_hit_count": len(hits),
        "kept_hit_count": len(filtered),
        "bait_count": len(rows),
        "total_unique_contig_count": len({hit.ctg_id for hit in deduped}),
        "bait_with_circular_contig_count": sum(1 for row in rows if int(row["circular_contig_count"]) > 0),
    }


def build_summary_params(args, actual_out: str) -> dict[str, object]:
    return {
        "command": "summarize",
        "status": "running",
        "started_at": now_iso(),
        "bait2contig_version": __version__,
        "python_version": ".".join(map(str, tuple(__import__("sys").version_info[:3]))),
        "out": str(Path(actual_out).resolve()),
        "gzip": bool(args.gzip),
        "hits": str(Path(args.hits).resolve()),
        "min_identity": args.min_identity if args.min_identity is not None else "NA",
        "min_coverage": args.min_coverage if args.min_coverage is not None else "NA",
        "min_aln_length": args.min_aln_length if args.min_aln_length is not None else "NA",
        "best_hit": bool(args.best_hit),
        "include_contigs": bool(args.include_contigs),
        "contig_sep": args.contig_sep,
    }


def summary_resume_params(start_params: dict[str, object]) -> dict[str, object]:
    keys = [
        "out",
        "gzip",
        "hits",
        "min_identity",
        "min_coverage",
        "min_aln_length",
        "best_hit",
        "include_contigs",
        "contig_sep",
    ]
    return {key: start_params[key] for key in keys}


def require_existing_file(path: str, option: str) -> None:
    file_path = Path(path)
    if not file_path.exists():
        raise SummaryError(f"{option} file not found: {path}")
    if not file_path.is_file():
        raise SummaryError(f"{option} must be a file: {path}")


def validate_output_path(path: Optional[str], option: str) -> None:
    if path is None:
        return
    output_path = Path(path)
    if output_path.exists() and output_path.is_dir():
        raise SummaryError(f"{option} must be a file path, not a directory: {path}")
    parent = output_path.parent
    if str(parent) and parent.exists() and not parent.is_dir():
        raise SummaryError(f"{option} parent path is not a directory: {parent}")


def validate_summary_args(args) -> None:
    if args.resume and args.rerun:
        raise SummaryError("--resume and --rerun cannot be used together.")
    require_existing_file(args.hits, "--hits")
    validate_output_path(args.out, "--out")
    validate_output_path(args.log, "--log")
    for label in ("min_identity", "min_coverage"):
        value = getattr(args, label)
        if value is not None and (value < 0 or value > 1):
            raise SummaryError(f"--{label.replace('_', '-')} must be between 0 and 1")
    if args.min_aln_length is not None and args.min_aln_length < 0:
        raise SummaryError("--min-aln-length must be at least 0")
    if args.monitor_interval < 1:
        raise SummaryError("--monitor-interval must be at least 1")


def check_existing_output(path: str, *, resume: bool, rerun: bool, force: bool) -> None:
    if not (resume or rerun or force) and Path(path).exists():
        raise SummaryError(
            f"output already exists: {path}\n"
            "Use --resume to reuse a completed run, --rerun to recompute, or --force to overwrite."
        )


def run_summarize(args) -> int:
    """Run bait2contig summarize from parsed argparse arguments."""

    validate_summary_args(args)
    actual_out = gzip_output_path(args.out, args.gzip)
    log_path = args.log or f"{actual_out}.log"
    start_params = build_summary_params(args, actual_out)

    if args.resume:
        with Logger(log_path, no_color=args.no_color, quiet=args.quiet, verbose=args.verbose) as logger:
            result = check_resume(
                log_path=log_path,
                command="summarize",
                output_path=actual_out,
                expected_params=summary_resume_params(start_params),
            )
            if result.ok:
                logger.info("Resume enabled: previous successful run found in log.")
                logger.info(f"Output verified: {Path(actual_out).name}")
                logger.info("Skipping summarize.")
                return 0
            logger.warn(f"Resume log found but {result.reason}")
            logger.info("Re-running summarize.")
    else:
        check_existing_output(actual_out, resume=args.resume, rerun=args.rerun, force=args.force)

    ensure_parent_dir(actual_out)
    logger = Logger(log_path, no_color=args.no_color, quiet=args.quiet, verbose=args.verbose)
    monitor = ResourceMonitor(args.monitor_interval, logger)
    monitor.start()
    start_time = time.time()
    try:
        logger.info(f"bait2contig version: {__version__}")
        logger.info("command: summarize")
        logger.info(f"output: {actual_out}")
        logger.info(f"log: {log_path}")
        logger.info(f"resume: {'enabled' if args.resume else 'disabled'}")
        logger.info(f"rerun: {'enabled' if args.rerun else 'disabled'}")
        logger.marker(START_MARKER, start_params)
        stats = summarize_file(
            hits_path=args.hits,
            out_path=actual_out,
            min_identity=args.min_identity,
            min_coverage=args.min_coverage,
            min_aln_length=args.min_aln_length,
            best_hit=args.best_hit,
            include_contigs=args.include_contigs,
            contig_sep=args.contig_sep,
        )
        resource_stats = monitor.stop()
        runtime = time.time() - start_time
        output_size = Path(actual_out).stat().st_size if Path(actual_out).exists() else 0
        output_lines = count_text_lines(actual_out) if Path(actual_out).exists() else 0
        logger.marker(
            DONE_MARKER,
            {
                "command": "summarize",
                "status": "success",
                "finished_at": now_iso(),
                "runtime_seconds": f"{runtime:.2f}",
                "exit_code": 0,
                "output": str(Path(actual_out).resolve()),
                "output_size": output_size,
                "output_lines": output_lines,
                "input_hit_count": stats["input_hit_count"],
                "kept_hit_count": stats["kept_hit_count"],
                "bait_count": stats["bait_count"],
                "total_unique_contig_count": stats["total_unique_contig_count"],
                "bait_with_circular_contig_count": stats["bait_with_circular_contig_count"],
                "peak_rss_mb": format_metric(resource_stats["peak_rss_mb"], digits=2),
                "mean_cpu_percent": format_metric(resource_stats["mean_cpu_percent"], digits=1),
                "max_cpu_percent": format_metric(resource_stats["max_cpu_percent"], digits=1),
            },
        )
        logger.done(f"wrote output: {actual_out}")
        logger.done("bait2contig summarize completed successfully")
        logger.done(f"runtime: {runtime:.2f} sec")
        logger.done(f"peak RSS: {format_metric(resource_stats['peak_rss_mb'], digits=2)} MB")
        logger.done(f"mean CPU: {format_metric(resource_stats['mean_cpu_percent'], digits=1)}%")
        logger.close()
        return 0
    except Exception as exc:
        resource_stats = monitor.stop()
        runtime = time.time() - start_time
        message = str(exc)
        logger.error(f"ERROR: {message}")
        logger.marker(
            FAILED_MARKER,
            {
                "command": "summarize",
                "status": "failed",
                "finished_at": now_iso(),
                "runtime_seconds": f"{runtime:.2f}",
                "exit_code": 1,
                "error_type": type(exc).__name__,
                "error": message,
                "peak_rss_mb": format_metric(resource_stats["peak_rss_mb"], digits=2),
                "mean_cpu_percent": format_metric(resource_stats["mean_cpu_percent"], digits=1),
                "max_cpu_percent": format_metric(resource_stats["max_cpu_percent"], digits=1),
            },
        )
        logger.close()
        raise LoggedSummaryError(message) from exc
