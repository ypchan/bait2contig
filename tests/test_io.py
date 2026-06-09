import pytest

from bait2contig.io import (
    gzip_output_path,
    infer_circular_from_header,
    open_text,
    parse_paf_line,
    read_fasta,
    read_lineage,
)


def test_open_text_plain(tmp_path):
    path = tmp_path / "plain.txt"
    with open_text(path, "wt") as handle:
        handle.write("bait2contig\n")
    with open_text(path, "rt") as handle:
        assert handle.read() == "bait2contig\n"


def test_open_text_gzip(tmp_path):
    path = tmp_path / "plain.txt.gz"
    with open_text(path, "wt") as handle:
        handle.write("bait2contig\n")
    with open_text(path, "rt") as handle:
        assert handle.read() == "bait2contig\n"


def test_gzip_output_path():
    assert gzip_output_path("bait2contig.hits.tsv", True) == "bait2contig.hits.tsv.gz"


def test_gzip_output_path_already_gz():
    assert gzip_output_path("bait2contig.hits.tsv.gz", True) == "bait2contig.hits.tsv.gz"


def test_fasta_header_id(tmp_path):
    path = tmp_path / "contigs.fa"
    path.write_text(">contig_001 len=15320 circular=true\nATGC\n", encoding="utf-8")
    records = read_fasta(path)
    assert set(records) == {"contig_001"}
    assert records["contig_001"].length == 4


def test_circular_header():
    assert infer_circular_from_header("contig_001 circular=true")
    assert infer_circular_from_header("contig_002 is_circular=true")
    assert infer_circular_from_header("contig_003 circular")
    assert infer_circular_from_header("contig_004 circ=true")
    assert not infer_circular_from_header("contig_005 circular=false")


def test_lineage_with_header(tmp_path):
    path = tmp_path / "lineage.tsv"
    path.write_text("bait_id\tlineage\nb1\td__A\nb2\td__B\n", encoding="utf-8")
    assert read_lineage(path) == {"b1": "d__A", "b2": "d__B"}


def test_lineage_without_header(tmp_path):
    path = tmp_path / "lineage.tsv"
    path.write_text("b1\td__A\nb2\td__B\n", encoding="utf-8")
    assert read_lineage(path) == {"b1": "d__A", "b2": "d__B"}


def test_lineage_duplicate_bait_id(tmp_path):
    path = tmp_path / "lineage.tsv"
    path.write_text("bait_id\tlineage\nb1\td__A\nb1\td__B\n", encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate bait_id"):
        read_lineage(path)


def test_paf_parse_basic():
    hit = parse_paf_line("bait1\t100\t5\t95\t+\tcontig1\t1000\t10\t100\t81\t90\t60")
    assert hit.bait_id == "bait1"
    assert hit.ctg_id == "contig1"
    assert hit.ctg_len == 1000
    assert hit.aln_length == 90
    assert hit.identity == pytest.approx(0.9)
    assert hit.cov_bait == pytest.approx(0.9)


def test_identity_and_cov_bait():
    hit = parse_paf_line("bait1\t200\t20\t180\t+\tcontig1\t1000\t10\t170\t144\t160\t60")
    assert hit.identity == pytest.approx(144 / 160)
    assert hit.cov_bait == pytest.approx(160 / 200)
