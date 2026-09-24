"""Property-based tests (hypothesis): descriptor invariants over random valid peptides."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
from hypothesis import given, settings, strategies as st
from pepx.descriptors import descriptor_vector, molecular_weight, AA_MASS, WATER
from pepx.alphabet import validate_sequence, encode_indices

PEPTIDE = st.text(alphabet="ACDEFGHIKLMNPQRSTVWY", min_size=1, max_size=60)

@given(PEPTIDE)
@settings(max_examples=300, deadline=None)
def test_mw_identity(seq):
    assert abs(molecular_weight(seq) - (sum(AA_MASS[a] for a in seq) - (len(seq)-1)*WATER)) < 1e-6

@given(PEPTIDE)
@settings(max_examples=300, deadline=None)
def test_descriptor_bounds(seq):
    d = descriptor_vector(seq)
    assert -4.6 <= d["gravy"] <= 4.6
    assert 0.0 <= d["aromaticity"] <= 1.0
    assert d["length"] == len(seq)
    assert d["mw"] > 0

@given(PEPTIDE)
@settings(max_examples=300, deadline=None)
def test_validate_encode_roundtrip(seq):
    assert validate_sequence(seq) == seq.upper()
    idx = encode_indices(seq)
    assert len(idx) == len(seq) and all(0 <= i < 20 for i in idx)
