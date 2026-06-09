from bait2contig.io import open_text, read_tsv
from bait2contig.summary import SummaryHit, summarize_file, summarize_hits


def hit(
    ctg_id,
    bait_id="bait1",
    identity=0.97,
    cov_bait=0.8,
    aln_length=100,
    ctg_len=1000,
    is_circular=False,
    lineage=None,
):
    return SummaryHit(
        ctg_id=ctg_id,
        bait_id=bait_id,
        identity=identity,
        cov_bait=cov_bait,
        aln_length=aln_length,
        ctg_len=ctg_len,
        is_circular=is_circular,
        lineage=lineage,
    )


def test_summarize_basic():
    rows = summarize_hits(
        [
            hit("c1", identity=0.99, cov_bait=0.9, aln_length=100, ctg_len=1000),
            hit("c2", identity=0.97, cov_bait=0.8, aln_length=80, ctg_len=2000, is_circular=True),
        ],
        include_lineage=False,
    )
    row = rows[0]
    assert row["contig_count"] == 2
    assert row["circular_contig_count"] == 1
    assert row["total_ctg_len"] == 3000
    assert row["mean_ctg_len"] == "1500.00"
    assert row["max_identity"] == "0.990000"


def test_summarize_deduplicate_same_bait_ctg():
    rows = summarize_hits(
        [
            hit("c1", identity=0.95, cov_bait=0.9, aln_length=100, ctg_len=1000),
            hit("c1", identity=0.99, cov_bait=0.8, aln_length=80, ctg_len=1000),
        ],
        include_lineage=False,
    )
    row = rows[0]
    assert row["contig_count"] == 1
    assert row["max_identity"] == "0.990000"


def test_summarize_with_lineage():
    rows = summarize_hits([hit("c1", lineage="d__A")], include_lineage=True)
    assert rows[0]["lineage"] == "d__A"


def test_summarize_best_hit():
    rows = summarize_hits(
        [hit("c1", identity=0.98), hit("c2", identity=0.99)],
        include_lineage=False,
        best_hit=True,
    )
    assert rows[0]["best_ctg_id"] == "c2"
    assert rows[0]["best_identity"] == "0.990000"


def test_summarize_include_contigs():
    rows = summarize_hits(
        [hit("c1", identity=0.98), hit("c2", identity=0.99), hit("c3", identity=0.97)],
        include_lineage=False,
        include_contigs=True,
    )
    assert rows[0]["contigs"] == "c2,c1,c3"


def test_summarize_empty_hits(tmp_path):
    hits = tmp_path / "hits.tsv"
    out = tmp_path / "summary.tsv"
    hits.write_text("ctg_id\tbait_id\tidentity\taln_length\tcov_bait\tctg_len\tis_circular\n", encoding="utf-8")
    summarize_file(hits_path=hits, out_path=out)
    columns, rows = read_tsv(out)
    assert rows == []
    assert columns[0] == "bait_id"
    assert "contig_count" in columns


def test_summarize_read_gzip_hits(tmp_path):
    hits = tmp_path / "hits.tsv.gz"
    out = tmp_path / "summary.tsv"
    with open_text(hits, "wt") as handle:
        handle.write("ctg_id\tbait_id\tidentity\taln_length\tcov_bait\tctg_len\tis_circular\n")
        handle.write("c1\tbait1\t0.990000\t100\t0.900000\t1000\tFalse\n")
    summarize_file(hits_path=hits, out_path=out, best_hit=True)
    _, rows = read_tsv(out)
    assert rows[0]["best_ctg_id"] == "c1"
