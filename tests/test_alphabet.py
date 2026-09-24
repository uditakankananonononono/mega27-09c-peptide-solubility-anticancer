import pytest

from pepx.alphabet import (InvalidSequenceError, Peptide, encode_indices,
                           validate_sequence)


def test_validate_uppercases_and_strips():
    assert validate_sequence("  klaklak ") == "KLAKLAK"


def test_validate_rejects_bad_residues():
    with pytest.raises(InvalidSequenceError):
        validate_sequence("ACDE1")
    with pytest.raises(InvalidSequenceError):
        validate_sequence("")


def test_validate_allows_x_when_requested():
    assert validate_sequence("ACXDE", allow_x=True) == "ACXDE"
    assert encode_indices("ACXDE", allow_x=True)[2] == 20


def test_peptide_immutable_validated():
    p = Peptide("klaklakklaklak", name="KLAKLAK")
    assert p.sequence == "KLAKLAKKLAKLAK"
    assert len(p) == 14
    with pytest.raises(InvalidSequenceError):
        Peptide("BJORK")
