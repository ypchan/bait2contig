from bait2contig.fasta_index import (
    FastaIndex,
    build_fasta_index,
    default_index_path,
    index_is_current,
    read_fasta_subset,
)


def test_build_fasta_index_and_fetch_records(tmp_path):
    fasta = tmp_path / "contigs.fa"
    fasta.write_text(
        ">c1 circular=true\nACGT\nACGT\n>c2\nTTTT\n",
        encoding="utf-8",
    )
    index = tmp_path / "contigs.fai"

    assert build_fasta_index(fasta, index, threads=2) == 2
    assert index_is_current(index, fasta)

    with FastaIndex(fasta, index) as fasta_index:
        info = fasta_index.get_info(["c1", "c2", "missing"])
        assert set(info) == {"c1", "c2"}
        assert info["c1"].length == 8
        assert info["c1"].is_circular
        assert info["c1"].line_bases == 4

        records = fasta_index.get_fasta_records(["c2", "c1"])
        assert records["c1"].sequence == "ACGTACGT"
        assert records["c2"].sequence == "TTTT"


def test_index_detects_changed_source(tmp_path):
    fasta = tmp_path / "contigs.fa"
    fasta.write_text(">c1\nACGT\n", encoding="utf-8")
    index = tmp_path / "contigs.fai"
    build_fasta_index(fasta, index)
    assert index_is_current(index, fasta)

    fasta.write_text(">c1\nACGTACGT\n", encoding="utf-8")
    assert not index_is_current(index, fasta)


def test_index_handles_last_record_without_final_newline(tmp_path):
    fasta = tmp_path / "contigs.fa"
    fasta.write_text(">c1\nAAAA\n>c2\nCCCC", encoding="utf-8")
    index = tmp_path / "contigs.fai"

    build_fasta_index(fasta, index, threads=2)
    with FastaIndex(fasta, index) as fasta_index:
        info = fasta_index.get_info(["c1", "c2"])
        assert info["c1"].length == 4
        assert info["c2"].length == 4
        records = fasta_index.get_fasta_records(["c2"])
        assert records["c2"].sequence == "CCCC"


def test_read_fasta_subset(tmp_path):
    fasta = tmp_path / "contigs.fa"
    fasta.write_text(">c1\nAAAA\n>c2 circular=true\nCCCC\n>c3\nGGGG\n", encoding="utf-8")

    records = read_fasta_subset(fasta, ["c2"])
    assert set(records) == {"c2"}
    assert records["c2"].sequence == "CCCC"
    assert records["c2"].is_circular


def test_default_index_path():
    assert default_index_path("contigs.fa") == "contigs.fa.bait2contig.fai"
