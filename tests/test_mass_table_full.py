"""Regression lock for the full 20-residue mass table (pyteomics cross-validation).
The modlAMP lock used KLAKLAKKLAKLAK (A/K/L only) and missed scrambled P,Q,V,W,Y."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
from pepx.descriptors import molecular_weight
from pyteomics import mass as pmass

def test_all_20_residue_masses_match_pyteomics():
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        assert abs(molecular_weight(aa) - pmass.calculate_mass(sequence=aa, average=True)) < 0.05, aa

def test_scrambled_residues_peptide():
    # spans every previously-wrong residue P,Q,V,W,Y
    seq = "PQVWYPQVWY"
    assert abs(molecular_weight(seq) - pmass.calculate_mass(sequence=seq, average=True)) < 0.05

def test_klaklak_lock_still_holds():
    assert abs(molecular_weight("KLAKLAKLAKLAK") - 1395.84) < 0.1
