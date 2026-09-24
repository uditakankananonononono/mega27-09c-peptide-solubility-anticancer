import pytest

from pepx.descriptors import (aliphatic_index, aromaticity, boman_index, gravy,
                              instability_index, isoelectric_point,
                              molecular_weight, net_charge)

# Published reference: KLA peptide (KLAKLAKKLAKLAK) is a known cationic ACP.
KLA = "KLAKLAKKLAKLAK"


def test_molecular_weight_water_subtraction():
    # A single residue: free amino acid mass; dipeptide: sum - 1 water.
    assert molecular_weight("A") == pytest.approx(89.09, abs=0.01)
    assert molecular_weight("AA") == pytest.approx(2 * 89.09 - 18.015, abs=0.01)


def test_net_charge_cationic_peptide():
    assert net_charge(KLA) > 5.0  # 6 K + N-term - C-term
    assert net_charge("EEEE") < -3.0


def test_pi_bisection_bounds_and_direction():
    assert isoelectric_point(KLA) > 9.5       # poly-lysine-rich: high pI
    assert isoelectric_point("EEEEEEEE") < 5  # poly-glutamate: low pI


def test_gravy_sign():
    assert gravy("IIIIII") > 3.0
    assert gravy("DDDDDD") < -3.0


def test_aromaticity():
    assert aromaticity("FWYAFWYA") == pytest.approx(0.75)


def test_aliphatic_index_ikai():
    # Pure alanine: X(A)=100 -> index 100.
    assert aliphatic_index("AAAAAAAA") == pytest.approx(100.0)


def test_boman_and_instability_ranges():
    assert -3.0 < boman_index(KLA) < 3.5
    assert instability_index(KLA) < 40  # stable peptide per Guruprasad cutoff


def test_instability_unstable_dipeptide():
    # P-P pair has DIWV 20.26; long poly-P should score unstable-ish vs poly-G.
    assert instability_index("PPPPPPPPPP") > instability_index("GGGGGGGGGG")
