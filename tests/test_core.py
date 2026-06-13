from pathlib import Path
from types import SimpleNamespace

import pytest

from bait2contig.cli import main
from bait2contig.core import (
    SearchError,
    SearchHit,
    distribute_minimap2_threads,
    extract_contigs,
    filter_hits,
    is_terminal_partial_hit,
    renamed_header,
    resolve_index_threads,
    resolve_minimap2_jobs,
    run_parallel_minimap2,
    select_best_per_bait,
    select_extraction_hits,
    split_bait_records,
)
from bait2contig.io import FastaRecord, open_text, read_fasta, read_tsv
from bait2contig.log import DONE_MARKER, FAILED_MARKER, START_MARKER, check_resume, parse_marker_blocks


class FakeResourceLogger:
    def __init__(self, verbose=False):
        self.verbose = verbose
        self.messages = []

    def resource(self, message):
        self.messages.append(message)


def hit(
    ctg_id,
    bait_id="bait1",
    identity=0.97,
    cov_bait=0.8,
    aln_length=100,
    ctg_len=1000,
    is_circular=False,
    lineage=None,
    bait_len=100,
    bait_start=0,
    bait_end=100,
    ctg_start=0,
    ctg_end=None,
):
    if ctg_end is None:
        ctg_end = ctg_len
    return SearchHit(
        ctg_id=ctg_id,
        bait_id=bait_id,
        identity=identity,
        cov_bait=cov_bait,
        aln_length=aln_length,
        ctg_len=ctg_len,
        is_circular=is_circular,
        bait_len=bait_len,
        bait_start=bait_start,
        bait_end=bait_end,
        ctg_start=ctg_start,
        ctg_end=ctg_end,
        lineage=lineage,
    )


def contigs():
    return {
        "c1": FastaRecord("c1", "c1 len=4", "ATGC", False),
        "c2": FastaRecord("c2", "c2 circular=true", "AAAA", True),
        "c3": FastaRecord("c3", "c3 len=4", "CCCC", False),
    }


def test_best_hit_ranking():
    hits = [
        hit("c1", identity=0.98, cov_bait=0.80, aln_length=100, ctg_len=1000),
        hit("c2", identity=0.98, cov_bait=0.90, aln_length=90, ctg_len=900),
        hit("c3", identity=0.97, cov_bait=1.00, aln_length=200, ctg_len=3000),
    ]
    assert select_best_per_bait(hits)[0].ctg_id == "c2"


def test_filter_hits():
    hits = [
        hit("c1", identity=0.96, cov_bait=0.9, aln_length=200),
        hit("c2", identity=0.99, cov_bait=0.7, aln_length=200),
        hit("c3", identity=0.99, cov_bait=0.9, aln_length=50),
        hit("c4", identity=0.99, cov_bait=0.9, aln_length=200),
    ]
    kept = filter_hits(hits, identity=0.97, coverage=0.8, min_aln_length=100)
    assert [item.ctg_id for item in kept] == ["c4"]


def test_terminal_filter_allows_full_coverage_middle_hit():
    item = hit("c1", cov_bait=1.0, bait_start=20, bait_end=120, bait_len=120, ctg_start=200, ctg_end=300)
    assert is_terminal_partial_hit(item, terminal_tolerance=5)


def test_terminal_filter_allows_partial_terminal_hit_with_tolerance():
    item = hit(
        "c1",
        cov_bait=0.95,
        bait_len=100,
        bait_start=5,
        bait_end=100,
        ctg_len=1000,
        ctg_start=900,
        ctg_end=1000,
    )
    assert is_terminal_partial_hit(item, terminal_tolerance=5)


def test_terminal_filter_rejects_partial_middle_hit():
    item = hit(
        "c1",
        cov_bait=0.60,
        bait_len=100,
        bait_start=20,
        bait_end=80,
        ctg_len=1000,
        ctg_start=200,
        ctg_end=800,
    )
    kept = filter_hits([item], identity=0.97, coverage=0.5, min_aln_length=0, terminal_tolerance=5)
    assert kept == []


def test_extract_contigs_all(tmp_path):
    out = tmp_path / "matched.fa"
    count = extract_contigs([hit("c1"), hit("c2")], contigs(), out, rename=False, include_lineage=False, dedup=True)
    assert count == 2
    records = read_fasta(out)
    assert set(records) == {"c1", "c2"}


def test_extract_contigs_best(tmp_path):
    selected = select_extraction_hits(
        [hit("c1", identity=0.98, cov_bait=0.8), hit("c2", identity=0.99, cov_bait=0.7)],
        mode="best",
        identity=0,
        coverage=0,
        min_aln_length=0,
    )
    out = tmp_path / "best.fa"
    count = extract_contigs(selected, contigs(), out, rename=False, include_lineage=False, dedup=True)
    assert count == 1
    assert set(read_fasta(out)) == {"c2"}


def test_extract_contigs_circular(tmp_path):
    selected = select_extraction_hits(
        [hit("c1", is_circular=False), hit("c2", is_circular=True)],
        mode="circular",
        identity=0,
        coverage=0,
        min_aln_length=0,
    )
    out = tmp_path / "circular.fa"
    extract_contigs(selected, contigs(), out, rename=False, include_lineage=False, dedup=True)
    assert set(read_fasta(out)) == {"c2"}


def test_extract_contigs_non_circular(tmp_path):
    selected = select_extraction_hits(
        [hit("c1", is_circular=False), hit("c2", is_circular=True)],
        mode="non-circular",
        identity=0,
        coverage=0,
        min_aln_length=0,
    )
    out = tmp_path / "non_circular.fa"
    extract_contigs(selected, contigs(), out, rename=False, include_lineage=False, dedup=True)
    assert set(read_fasta(out)) == {"c1"}


def test_extract_rename(tmp_path):
    item = hit("c1", bait_id="baitA", identity=0.991234, cov_bait=0.982394, aln_length=1450, ctg_len=32781)
    header = renamed_header(item)
    assert "bait=baitA" in header
    assert "identity=0.991234" in header
    assert "cov_bait=0.982394" in header
    assert "aln_length=1450" in header
    assert "ctg_len=32781" in header
    assert "circular=False" in header


def test_extract_include_lineage_requires_rename(tmp_path):
    with pytest.raises(SearchError, match="requires --extract-rename"):
        extract_contigs([hit("c1")], contigs(), tmp_path / "x.fa", rename=False, include_lineage=True, dedup=True)


class FakeLogger:
    def __init__(self):
        self.messages = []
        self.verbose = False

    def info(self, message):
        self.messages.append(message)

    def warn(self, message):
        self.messages.append(message)


def test_resource_monitor_periodic_logs_require_verbose():
    from bait2contig.log import ResourceMonitor

    logger = FakeResourceLogger(verbose=False)
    monitor = ResourceMonitor(1, logger)
    monitor.sample(write_log=True)
    assert logger.messages == []

    verbose_logger = FakeResourceLogger(verbose=True)
    verbose_monitor = ResourceMonitor(1, verbose_logger)
    verbose_monitor.sample(write_log=True)
    assert verbose_logger.messages


def test_logger_uses_short_human_timestamp(tmp_path):
    from bait2contig.log import Logger

    log = tmp_path / "run.log"
    logger = Logger(log, no_color=True, quiet=True)
    logger.info("hello")
    logger.close()

    line = log.read_text(encoding="utf-8").splitlines()[0]
    assert line[4] == "-"
    assert line[10] == " "
    assert "T" not in line.split()[0]
    assert "+08:00" not in line


def test_no_extract_dedup_warning(tmp_path):
    logger = FakeLogger()
    extract_contigs(
        [hit("c1"), hit("c1")],
        contigs(),
        tmp_path / "x.fa",
        rename=False,
        include_lineage=False,
        dedup=False,
        logger=logger,
    )
    assert any("--no-extract-dedup" in message for message in logger.messages)


def test_auto_index_threads_follow_mapping_threads():
    args = SimpleNamespace(index_threads=0, threads=32)

    assert resolve_index_threads(args) == 32


def test_explicit_index_threads_override_mapping_threads():
    args = SimpleNamespace(index_threads=4, threads=32)

    assert resolve_index_threads(args) == 4


def test_resolve_minimap2_jobs_caps_to_baits_and_threads():
    args = SimpleNamespace(minimap2_jobs=64, threads=32)

    assert resolve_minimap2_jobs(args, bait_count=10) == 10
    assert resolve_minimap2_jobs(args, bait_count=64) == 32


def test_distribute_minimap2_threads_uses_total_budget():
    assert distribute_minimap2_threads(32, 3) == [11, 11, 10]
    assert distribute_minimap2_threads(2, 8) == [1, 1]


def test_split_bait_records_keeps_contiguous_chunks():
    records = [
        FastaRecord("b1", "b1", "A"),
        FastaRecord("b2", "b2", "C"),
        FastaRecord("b3", "b3", "G"),
        FastaRecord("b4", "b4", "T"),
    ]

    chunks = split_bait_records(records, 2)

    assert [[record.id for record in chunk] for chunk in chunks] == [["b1", "b2"], ["b3", "b4"]]


def test_run_parallel_minimap2_combines_chunks_in_order(tmp_path, monkeypatch):
    records = [
        FastaRecord("b1", "b1", "A"),
        FastaRecord("b2", "b2", "C"),
        FastaRecord("b3", "b3", "G"),
        FastaRecord("b4", "b4", "T"),
    ]
    tmp_paf = tmp_path / "combined.paf"

    def fake_execute(command, paf_path, monitor):
        bait_path = Path(command[-2])
        headers = [
            line[1:].strip()
            for line in bait_path.read_text(encoding="utf-8").splitlines()
            if line.startswith(">")
        ]
        Path(paf_path).write_text("".join(f"{header}\tthreads={command[4]}\n" for header in headers), encoding="utf-8")
        return 0, ""

    monkeypatch.setattr("bait2contig.core.execute_minimap2_command", fake_execute)

    run_parallel_minimap2(
        executable="minimap2",
        preset="asm10",
        total_threads=4,
        jobs=2,
        contigs="contigs.fa",
        bait_records=records,
        tmp_paf=str(tmp_paf),
        tmp_dir=tmp_path,
        logger=FakeLogger(),
        monitor=SimpleNamespace(set_child_pid=lambda pid: None),
    )

    assert tmp_paf.read_text(encoding="utf-8").splitlines() == [
        "b1\tthreads=2",
        "b2\tthreads=2",
        "b3\tthreads=2",
        "b4\tthreads=2",
    ]


def test_search_maps_contigs_as_query_and_adds_lineage(tmp_path):
    contigs = tmp_path / "contigs.fa"
    bait = tmp_path / "bait.fa"
    lineage = tmp_path / "lineage.tsv"
    out = tmp_path / "hits.tsv"
    fake_minimap2 = tmp_path / "minimap2"

    contigs.write_text(">ctg1\n" + "A" * 1000 + "\n", encoding="utf-8")
    bait.write_text(">bait1\n" + "A" * 1000 + "\n", encoding="utf-8")
    lineage.write_text("bait_id\tlineage\nbait1\td__A;p__B\n", encoding="utf-8")
    fake_minimap2.write_text(
        "\n".join(
            [
                "#!/usr/bin/env python3",
                "import pathlib",
                "import sys",
                "if '--version' in sys.argv:",
                "    print('fake-minimap2')",
                "    raise SystemExit(0)",
                "target = pathlib.Path(sys.argv[-2]).name",
                "query = pathlib.Path(sys.argv[-1]).name",
                "if target != 'bait.fa' or query != 'contigs.fa':",
                "    print(f'unexpected order: target={target} query={query}', file=sys.stderr)",
                "    raise SystemExit(2)",
                "print('ctg1\\t1000\\t10\\t900\\t+\\tbait1\\t1000\\t0\\t990\\t980\\t990\\t60')",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    fake_minimap2.chmod(0o755)

    exit_code = main(
        [
            "search",
            "--contigs",
            str(contigs),
            "--bait",
            str(bait),
            "--lineage",
            str(lineage),
            "--out",
            str(out),
            "--identity",
            "0.9",
            "--coverage",
            "0.9",
            "--no-terminal-filter",
            "--minimap2",
            str(fake_minimap2),
            "--quiet",
        ]
    )

    assert exit_code == 0
    columns, rows = read_tsv(out)
    assert "lineage" in columns
    assert rows == [
        {
            "ctg_id": "ctg1",
            "bait_id": "bait1",
            "identity": "0.989899",
            "aln_length": "990",
            "cov_bait": "0.990000",
            "ctg_len": "1000",
            "is_circular": "False",
            "lineage": "d__A;p__B",
        }
    ]


def write_resume_log(path, command, output, params, marker=DONE_MARKER):
    lines = [START_MARKER, f"command={command}", "status=running"]
    lines.extend(f"{key}={value}" for key, value in params.items())
    lines.extend(
        [
            marker,
            f"command={command}",
            "status=success" if marker == DONE_MARKER else "status=failed",
            "exit_code=0" if marker == DONE_MARKER else "exit_code=1",
            f"output={output}",
            "output_size=12",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_log_parser_done_success(tmp_path):
    out = tmp_path / "bait2contig.hits.tsv"
    out.write_text("ctg_id\n", encoding="utf-8")
    log = tmp_path / "run.log"
    params = {"out": str(out.resolve()), "gzip": "false", "identity": "0.97"}
    write_resume_log(log, "search", str(out.resolve()), params)
    blocks = parse_marker_blocks(log)
    assert blocks[-1].marker == DONE_MARKER
    result = check_resume(log_path=log, command="search", output_path=out, expected_params=params)
    assert result.ok


def test_log_parser_incomplete_start_only(tmp_path):
    out = tmp_path / "bait2contig.hits.tsv"
    out.write_text("ctg_id\n", encoding="utf-8")
    log = tmp_path / "run.log"
    log.write_text(f"{START_MARKER}\ncommand=search\nout={out.resolve()}\n", encoding="utf-8")
    result = check_resume(log_path=log, command="search", output_path=out, expected_params={"out": str(out.resolve())})
    assert not result.ok
    assert "incomplete" in result.reason


def test_log_parser_failed(tmp_path):
    out = tmp_path / "bait2contig.hits.tsv"
    out.write_text("ctg_id\n", encoding="utf-8")
    log = tmp_path / "run.log"
    params = {"out": str(out.resolve()), "gzip": "false"}
    write_resume_log(log, "search", str(out.resolve()), params, marker=FAILED_MARKER)
    result = check_resume(log_path=log, command="search", output_path=out, expected_params=params)
    assert not result.ok
    assert "failed" in result.reason


def test_resume_rejects_changed_identity(tmp_path):
    out = tmp_path / "bait2contig.hits.tsv"
    out.write_text("ctg_id\n", encoding="utf-8")
    log = tmp_path / "run.log"
    params = {"out": str(out.resolve()), "gzip": "false", "identity": "0.97"}
    write_resume_log(log, "search", str(out.resolve()), params)
    result = check_resume(
        log_path=log,
        command="search",
        output_path=out,
        expected_params={"out": str(out.resolve()), "gzip": "false", "identity": "0.99"},
    )
    assert not result.ok
    assert "identity 0.97 -> 0.99" in result.reason


def test_resume_accepts_matching_successful_log(tmp_path):
    out = tmp_path / "bait2contig.hits.tsv"
    out.write_text("ctg_id\n", encoding="utf-8")
    log = tmp_path / "run.log"
    params = {"out": str(out.resolve()), "gzip": "false", "identity": "0.97"}
    write_resume_log(log, "search", str(out.resolve()), params)
    assert check_resume(log_path=log, command="search", output_path=out, expected_params=params).ok


def test_resume_with_gzip_output(tmp_path):
    out = tmp_path / "bait2contig.hits.tsv.gz"
    with open_text(out, "wt") as handle:
        handle.write("ctg_id\n")
    log = tmp_path / "run.log"
    params = {"out": str(out.resolve()), "gzip": "true"}
    write_resume_log(log, "search", str(out.resolve()), params)
    assert check_resume(log_path=log, command="search", output_path=out, expected_params=params).ok


def test_cli_help(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--help"])
    assert exc.value.code == 0
    assert "search" in capsys.readouterr().out


def test_cli_search_help(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["search", "--help"])
    assert exc.value.code == 0
    output = capsys.readouterr().out
    assert "Required arguments" in output
    assert "Contig extraction arguments" in output
    assert "--terminal-tolerance" in output
    assert "Minimum alignment identity. (default: 0.97)" in output
    assert "Total minimap2 thread budget. (default: 8)" in output
    assert "Parallel minimap2 processes for splitting bait FASTA" in output
    assert "--bait-index is not used. (default: 1)" in output
    assert "--bait-index" in output
    assert "Text FASTA index path. Default: <contigs>.bait2contig.fai." in output
    assert "Threads for building plain FASTA indexes. Use 0 to follow --threads." in output
    assert "Minimum identity for extracted contigs. Default: --identity." in output
    assert "Input contig FASTA. (default: None)" not in output


def test_cli_summarize_help(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["summarize", "--help"])
    assert exc.value.code == 0
    output = capsys.readouterr().out
    assert "Summary arguments" in output
    assert "Resume and output arguments" in output
    assert "Separator for contig lists. (default: ,)" in output
    assert "Additional identity filter. Default: no additional filter." in output


def test_cli_missing_command(capsys):
    assert main([]) == 2
    error = capsys.readouterr().err
    assert "missing command" in error
    assert "search, summarize" in error


def test_cli_unknown_command_suggestion(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["seach"])
    assert exc.value.code == 2
    error = capsys.readouterr().err
    assert "invalid value for <command>" in error
    assert "Did you mean 'search'?" in error


def test_cli_unknown_option_suggestion(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["search", "--contigs", "contigs.fa", "--bait", "bait.fa", "--out", "out.tsv", "--identit", "0.9"])
    assert exc.value.code == 2
    error = capsys.readouterr().err
    assert "unrecognized argument" in error
    assert "Did you mean --identity?" in error


def test_cli_missing_required_arguments_has_example(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["search"])
    assert exc.value.code == 2
    error = capsys.readouterr().err
    assert "missing required argument" in error
    assert "bait2contig search --contigs contigs.fa --bait bait.fa --out bait2contig.hits.tsv" in error


def test_search_missing_input_file_is_reported(capsys, tmp_path):
    bait = tmp_path / "bait.fa"
    bait.write_text(">bait1\nATGC\n", encoding="utf-8")
    exit_code = main(["search", "--contigs", str(tmp_path / "missing.fa"), "--bait", str(bait), "--out", str(tmp_path / "out.tsv")])
    assert exit_code == 1
    assert "--contigs file not found" in capsys.readouterr().err


def test_summarize_missing_hits_file_is_reported(capsys, tmp_path):
    exit_code = main(["summarize", "--hits", str(tmp_path / "missing.tsv"), "--out", str(tmp_path / "summary.tsv")])
    assert exit_code == 1
    assert "--hits file not found" in capsys.readouterr().err
