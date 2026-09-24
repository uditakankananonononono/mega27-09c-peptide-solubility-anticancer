import pytest

from pepx.fasta import read_fasta, write_fasta


def test_roundtrip(tmp_path):
    p = tmp_path / "x.fa"
    n = write_fasta(iter([("h1 desc", "ACDE"), ("h2", "FGHIK")]), p)
    assert n == 2
    recs = list(read_fasta(p))
    assert recs == [("h1 desc", "ACDE"), ("h2", "FGHIK")]


def test_wrapping(tmp_path):
    p = tmp_path / "y.fa"
    write_fasta(iter([("h", "A" * 100)]), p, wrap=30)
    header, seq = list(read_fasta(p))[0]
    assert seq == "A" * 100


def test_malformed_raises(tmp_path):
    p = tmp_path / "z.fa"
    p.write_text("ACDE\n>h\nFG\n")
    with pytest.raises(ValueError):
        list(read_fasta(p))
