import pytest

from bait2contig.cli import main
from bait2contig.core import (
    SearchError,
    SearchHit,
    extract_contigs,
    filter_hits,
    is_terminal_partial_hit,
    renamed_header,
    select_best_per_bait,
    select_extraction_hits,
)
from bait2contig.io import FastaRecord, open_text, read_fasta
from bait2contig.log import DONE_MARKER, FAILED_MARKER, START_MARKER, check_resume, parse_marker_blocks


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

    def warn(self, message):
        self.messages.append(message)


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


def test_cli_summarize_help(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["summarize", "--help"])
    assert exc.value.code == 0
    output = capsys.readouterr().out
    assert "Summary arguments" in output
    assert "Resume and output arguments" in output


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
