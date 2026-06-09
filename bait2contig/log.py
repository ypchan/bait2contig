"""Logging, resource monitoring, and resume parsing for bait2contig."""

from __future__ import annotations

import os
import sys
import threading
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from .io import ensure_parent_dir, open_text

try:  # Optional dependency.
    import psutil  # type: ignore
except Exception:  # pragma: no cover - depends on the runtime environment.
    psutil = None  # type: ignore

try:
    import resource
except Exception:  # pragma: no cover - Windows fallback.
    resource = None  # type: ignore


START_MARKER = "[BAIT2CONTIG_START]"
DONE_MARKER = "[BAIT2CONTIG_DONE]"
FAILED_MARKER = "[BAIT2CONTIG_FAILED]"
STANDARD_MARKERS = {START_MARKER, DONE_MARKER, FAILED_MARKER}


@dataclass(frozen=True)
class ResumeResult:
    ok: bool
    reason: str


@dataclass(frozen=True)
class MarkerBlock:
    marker: str
    values: Dict[str, str]


class Logger:
    """Write grep-friendly logs to stderr and a plain-text log file."""

    COLORS = {
        "INFO": "\033[36m",
        "WARN": "\033[33m",
        "ERROR": "\033[31m",
        "DONE": "\033[32m",
        "RESOURCE": "\033[2m",
    }
    RESET = "\033[0m"

    def __init__(
        self,
        log_path: str | Path,
        *,
        no_color: bool = False,
        quiet: bool = False,
        verbose: bool = False,
        append: bool = True,
    ) -> None:
        self.log_path = str(log_path)
        ensure_parent_dir(self.log_path)
        self._handle = open(self.log_path, "a" if append else "w", encoding="utf-8", newline="")
        self.quiet = quiet
        self.verbose = verbose
        self.color = (not no_color) and sys.stderr.isatty()

    def close(self) -> None:
        if not self._handle.closed:
            self._handle.close()

    def _screen_enabled(self, level: str) -> bool:
        if self.quiet and level not in {"WARN", "ERROR"}:
            return False
        return True

    def emit(self, level: str, message: str) -> None:
        tag = f"[{level}]".ljust(11)
        plain_line = f"{tag}{message}"
        self._handle.write(f"{plain_line}\n")
        self._handle.flush()
        if self._screen_enabled(level):
            if self.color and level in self.COLORS:
                sys.stderr.write(f"{self.COLORS[level]}{plain_line}{self.RESET}\n")
            else:
                sys.stderr.write(f"{plain_line}\n")
            sys.stderr.flush()

    def info(self, message: str) -> None:
        self.emit("INFO", message)

    def warn(self, message: str) -> None:
        self.emit("WARN", message)

    def error(self, message: str) -> None:
        self.emit("ERROR", message)

    def done(self, message: str) -> None:
        self.emit("DONE", message)

    def resource(self, message: str) -> None:
        self.emit("RESOURCE", message)

    def marker(self, marker: str, values: Dict[str, object]) -> None:
        self._handle.write(f"{marker}\n")
        for key, value in values.items():
            self._handle.write(f"{key}={format_log_value(value)}\n")
        self._handle.flush()

    def __enter__(self) -> "Logger":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()


class ResourceMonitor:
    """Collect process resource usage without making monitoring mandatory."""

    def __init__(self, interval: int, logger: Logger) -> None:
        self.interval = max(1, int(interval))
        self.logger = logger
        self.start_time = time.time()
        self._stop = threading.Event()
        self._thread: Optional[threading.Thread] = None
        self._child_pid: Optional[int] = None
        self._cpu_samples: List[float] = []
        self._peak_rss_mb: Optional[float] = None
        self._process = psutil.Process(os.getpid()) if psutil is not None else None
        if self._process is not None:
            try:
                self._process.cpu_percent(None)
            except Exception:
                pass

    def set_child_pid(self, pid: int) -> None:
        self._child_pid = pid
        if psutil is not None:
            try:
                psutil.Process(pid).cpu_percent(None)
            except Exception:
                pass

    def start(self) -> None:
        self._thread = threading.Thread(target=self._run, name="bait2contig-resource-monitor", daemon=True)
        self._thread.start()

    def stop(self) -> Dict[str, Optional[float]]:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=1)
        self.sample(write_log=False)
        return {
            "peak_rss_mb": self._peak_rss_mb,
            "mean_cpu_percent": mean(self._cpu_samples),
            "max_cpu_percent": max(self._cpu_samples) if self._cpu_samples else None,
        }

    def _run(self) -> None:
        while not self._stop.wait(self.interval):
            self.sample(write_log=True)

    def sample(self, *, write_log: bool) -> None:
        try:
            rss_mb, cpu_percent = self._collect()
            if rss_mb is not None:
                if self._peak_rss_mb is None or rss_mb > self._peak_rss_mb:
                    self._peak_rss_mb = rss_mb
            if cpu_percent is not None:
                self._cpu_samples.append(cpu_percent)
            if write_log:
                elapsed = time.time() - self.start_time
                self.logger.resource(
                    "elapsed={:.1f}s rss_mb={} cpu_percent={}".format(
                        elapsed,
                        format_metric(rss_mb, digits=1),
                        format_metric(cpu_percent, digits=1),
                    )
                )
        except Exception:
            return

    def _collect(self) -> tuple[Optional[float], Optional[float]]:
        if psutil is not None and self._process is not None:
            rss = 0
            cpu = 0.0
            processes = [self._process]
            try:
                processes.extend(self._process.children(recursive=True))
            except Exception:
                pass
            if self._child_pid is not None:
                try:
                    child = psutil.Process(self._child_pid)
                    if child not in processes:
                        processes.append(child)
                except Exception:
                    pass
            for proc in processes:
                try:
                    rss += proc.memory_info().rss
                    cpu += proc.cpu_percent(None)
                except Exception:
                    continue
            return (rss / (1024 * 1024), cpu)
        return (resource_peak_rss_mb(), None)


def now_iso() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def format_log_value(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def format_metric(value: Optional[float], digits: int = 2) -> str:
    if value is None:
        return "NA"
    return f"{value:.{digits}f}"


def mean(values: Iterable[float]) -> Optional[float]:
    values = list(values)
    if not values:
        return None
    return sum(values) / len(values)


def resource_peak_rss_mb() -> Optional[float]:
    if resource is None:
        return None
    try:
        own = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        raw = max(own, child)
        if sys.platform == "darwin":
            return raw / (1024 * 1024)
        return raw / 1024
    except Exception:
        return None


def parse_marker_blocks(log_path: str | Path) -> List[MarkerBlock]:
    """Parse standardized marker blocks from a plain-text log."""

    blocks: List[MarkerBlock] = []
    current: Optional[MarkerBlock] = None
    with open_text(log_path, "rt") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line in STANDARD_MARKERS:
                current = MarkerBlock(marker=line, values={})
                blocks.append(current)
                continue
            if current is not None and "=" in line and not line.startswith("["):
                key, value = line.split("=", 1)
                current.values[key] = value
    return blocks


def check_resume(
    *,
    log_path: str | Path,
    command: str,
    output_path: str | Path,
    expected_params: Dict[str, object],
) -> ResumeResult:
    """Return whether a previous successful run can be resumed."""

    log_path = str(log_path)
    output_path = str(Path(output_path).resolve())
    output_file = Path(output_path)
    if not Path(log_path).exists():
        return ResumeResult(False, "resume log not found")
    if not output_file.exists() or output_file.stat().st_size <= 0:
        return ResumeResult(False, "output missing or empty")
    try:
        blocks = parse_marker_blocks(log_path)
    except Exception as exc:
        return ResumeResult(False, f"resume log is malformed: {exc}")
    if not blocks:
        return ResumeResult(False, "resume log has no standardized marker blocks")

    latest_index: Optional[int] = None
    for index in range(len(blocks) - 1, -1, -1):
        if blocks[index].values.get("command") == command:
            latest_index = index
            break
    if latest_index is None:
        return ResumeResult(False, f"no previous {command} marker found")

    latest = blocks[latest_index]
    if latest.marker == START_MARKER:
        return ResumeResult(False, "previous run is incomplete")
    if latest.marker == FAILED_MARKER:
        return ResumeResult(False, "previous run failed")
    if latest.marker != DONE_MARKER:
        return ResumeResult(False, "latest marker is not a successful run")
    if latest.values.get("status") != "success" or latest.values.get("exit_code") != "0":
        return ResumeResult(False, "previous run did not finish successfully")
    if str(Path(latest.values.get("output", "")).resolve()) != output_path:
        return ResumeResult(False, "previous output path does not match")
    try:
        if int(latest.values.get("output_size", "0")) <= 0:
            return ResumeResult(False, "previous output size is zero")
    except ValueError:
        return ResumeResult(False, "previous output size is malformed")

    start_block: Optional[MarkerBlock] = None
    for index in range(latest_index - 1, -1, -1):
        block = blocks[index]
        if block.marker == START_MARKER and block.values.get("command") == command:
            start_block = block
            break
    if start_block is None:
        return ResumeResult(False, "matching start marker was not found")

    for key, expected in expected_params.items():
        expected_value = format_log_value(expected)
        observed = start_block.values.get(key)
        if observed is None:
            return ResumeResult(False, f"resume log missing parameter: {key}")
        if observed != expected_value:
            return ResumeResult(False, f"parameters changed: {key} {observed} -> {expected_value}")

    return ResumeResult(True, "previous successful run found in log")
