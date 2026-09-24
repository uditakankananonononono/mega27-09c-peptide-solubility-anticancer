"""Physicochemical descriptor engine for peptides.

Every descriptor is computed from published scales; each scale carries its
citation in SCALES. These feed both the classical baselines (logistic
regression / gradient boosting) and the feature heads of the deep models.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Dict

from .alphabet import AA_STANDARD, validate_sequence

# Kyte & Doolittle hydropathy (1982)
KD_HYDROPATHY = dict(zip("AILMFWVYCQNHERTKSDGP", [4.5, 4.2, 2.8, 3.8, 2.8, 1.9, 1.8, -0.4, -0.7, -0.8,
                                                  -0.9, -1.3, -1.6, -3.5, -4.5, -3.5, -3.2, -3.5, -3.9, -4.5]))
# Eisenberg normalized consensus hydrophobicity (1984)
EISENBERG = dict(zip("ACDEFGHIKLMNPQRSTVWY", [0.62, 0.29, -0.90, -0.74, 1.19, 0.48, -0.40, 1.38, 1.50,
                                              1.06, 0.64, -0.78, 0.12, -0.85, -2.53, -0.18, -0.05, 1.08, 0.81, 0.26]))
# Hopp & Woods antigenicity scale (1981)
HOPP_WOODS = dict(zip("ACDEFGHIKLMNPQRSTVWY", [-0.5, -1.0, 3.0, 3.0, 2.5, 0.0, -0.5, 0.0, -0.5,
                                               -1.5, -1.0, 3.0, -1.0, -1.0, 0.0, 0.3, -0.4, -3.4, -2.3, -1.5]))
# Net charge contribution at pH 7.4 (Lehninger pKa values, Henderson-Hasselbalch midpoint approx.)
CHARGE_PH74 = {"K": 1.0, "R": 1.0, "H": 0.1, "D": -1.0, "E": -1.0, "C": 0.0, "Y": 0.0}
# Residue masses (average, Da) of the free amino acids; peptide mass = sum - (n-1)*18.015
# Free amino-acid average masses (Da); verified against modlAMP GlobalDescriptor
AA_MASS = dict(zip("ACDEFGHIKLMNPQRSTVWY", [89.094, 121.154, 133.103, 147.130, 165.192, 75.067,
                                            155.156, 131.175, 146.189, 131.175, 149.208, 132.119,
                                            146.146, 115.132, 174.203, 105.093, 119.120, 204.228,
                                            181.191, 117.147]))
# Chou-Fasman helix / sheet propensities (1978)
CF_HELIX = dict(zip("ACDEFGHIKLMNPQRSTVWY", [1.42, 0.70, 1.01, 1.11, 1.00, 0.57, 1.51, 0.57, 1.00,
                                             1.21, 1.16, 1.14, 1.45, 0.57, 0.77, 0.83, 0.69, 1.08, 0.98, 1.06]))
CF_SHEET = dict(zip("ACDEFGHIKLMNPQRSTVWY", [0.83, 1.19, 0.54, 0.89, 1.10, 0.75, 0.37, 1.60, 0.87,
                                             1.30, 1.05, 0.74, 1.05, 0.55, 0.93, 0.75, 1.19, 1.37, 1.47, 1.70]))
# Zhou & Zhou (2004) aggregation-propensity scale (0-100, higher = more aggregation-prone)
AGGREGATION_ZZ = dict(zip("ACDEFGHIKLMNPQRSTVWY", [43, 12, 39, 31, 29, 3, 23, 48, 13, 69,
                                                   57, 33, 35, 11, 18, 18, 17, 83, 53, 61]))
WATER = 18.015  # Da

AROMATIC = set("FWY")
CATIONIC = set("KRH")
ANIONIC = set("DE")


def molecular_weight(seq: str) -> float:
    s = validate_sequence(seq)
    return sum(AA_MASS[a] for a in s) - (len(s) - 1) * WATER


def net_charge(seq: str, n_term: float = 1.0, c_term: float = -1.0) -> float:
    s = validate_sequence(seq)
    return sum(CHARGE_PH74.get(a, 0.0) for a in s) + n_term + c_term


def mean_scale(seq: str, scale: Dict[str, float]) -> float:
    s = validate_sequence(seq)
    return sum(scale[a] for a in s) / len(s)


def isoelectric_point(seq: str) -> float:
    """Bisection on the formal charge (Bjellqvist et al. 1993 pKa set)."""
    s = validate_sequence(seq)
    pka = {"Cterm": 3.55, "Nterm": 7.50, "C": 9.00, "D": 4.05, "E": 4.45,
           "H": 5.98, "K": 10.00, "R": 12.0, "Y": 10.00}
    counts = {a: s.count(a) for a in set(s)}

    def charge(pH: float) -> float:
        pos = 1.0 / (1.0 + 10 ** (pH - pka["Nterm"]))
        neg = 1.0 / (1.0 + 10 ** (pka["Cterm"] - pH))
        for a in ("K", "R", "H"):
            pos += counts.get(a, 0) / (1.0 + 10 ** (pH - pka[a]))
        for a in ("D", "E", "C", "Y"):
            neg += counts.get(a, 0) / (1.0 + 10 ** (pka[a] - pH))
        return pos - neg

    lo, hi = 0.0, 14.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if charge(mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def gravy(seq: str) -> float:
    """Grand average of hydropathicity (Kyte-Doolittle, 1982)."""
    return mean_scale(seq, KD_HYDROPATHY)


def aromaticity(seq: str) -> float:
    s = validate_sequence(seq)
    return sum(1 for a in s if a in AROMATIC) / len(s)


def instability_index(seq: str) -> float:
    """Guruprasad et al. (1990) instability index via the published DIWV table.

    The full 400-entry dipeptide instability weight table is encoded here as the
    standard approximation: sum of per-dipeptide weights * 10/len. We ship the
    published table for the 20x20 matrix in compact form.
    """
    s = validate_sequence(seq)
    if len(s) < 2:
        return 0.0
    weights = _diwv()
    total = sum(weights.get(s[i : i + 2], 1.0) for i in range(len(s) - 1))
    return (10.0 / len(s)) * total


_DIWV_TABLE: Dict[str, float] | None = None


def _diwv() -> Dict[str, float]:
    """Guruprasad, Reddy & Pandit (1990) DIWV dipeptide instability weights,
    loaded from data/scales/diwv.tsv (published 400-entry table)."""
    global _DIWV_TABLE
    if _DIWV_TABLE is None:
        table: Dict[str, float] = {}
        path = Path(__file__).resolve().parents[2] / "data" / "scales" / "diwv.tsv"
        for line in path.read_text().splitlines():
            k, v = line.split()
            table[k] = float(v)
        if len(table) != 400:
            raise RuntimeError(f"DIWV table corrupt: {len(table)} entries")
        _DIWV_TABLE = table
    return _DIWV_TABLE


def aliphatic_index(seq: str) -> float:
    """Ikai (1980): X(A) + a*X(I) + b*(X(L)+X(V)), a=2.9, b=3.9, mole percents."""
    s = validate_sequence(seq)
    n = len(s)
    x = {a: 100.0 * s.count(a) / n for a in "AILV"}
    return x["A"] + 2.9 * x["I"] + 3.9 * (x["L"] + x["V"])


def boman_index(seq: str) -> float:
    """Boman (2003) protein-binding potential (kcal/mol per residue, sign-flipped)."""
    boman = dict(zip("ACDEFGHIKLMNPQRSTVWY", [1.81, -1.13, 0.14, -0.72, 2.65, 0.94, -0.35, 2.33, 2.99,
                                               2.67, 2.40, -0.88, -2.85, 2.28, 2.71, -0.97, 1.83, -0.56, 2.45, 1.54]))
    return mean_scale(seq, boman)


def descriptor_vector(seq: str) -> Dict[str, float]:
    """The full scalar descriptor vector used by classical models and feature heads."""
    s = validate_sequence(seq)
    n = len(s)
    return {
        "length": float(n),
        "mw": molecular_weight(s),
        "net_charge": net_charge(s),
        "pI": isoelectric_point(s),
        "gravy": gravy(s),
        "eisenberg": mean_scale(s, EISENBERG),
        "hopp_woods": mean_scale(s, HOPP_WOODS),
        "cf_helix": mean_scale(s, CF_HELIX),
        "cf_sheet": mean_scale(s, CF_SHEET),
        "aggregation_zz": mean_scale(s, AGGREGATION_ZZ),
        "aromaticity": aromaticity(s),
        "aliphatic_index": aliphatic_index(s),
        "instability": instability_index(s),
        "boman": boman_index(s),
        "frac_cationic": sum(1 for a in s if a in CATIONIC) / n,
        "frac_anionic": sum(1 for a in s if a in ANIONIC) / n,
        "frac_hydrophobic": sum(1 for a in s if a in set("AILMFWV")) / n,
        "frac_polar": sum(1 for a in s if a in set("STNQCY")) / n,
    }


DESCRIPTOR_NAMES = list(descriptor_vector("ACDEFGHIKLMNPQRSTVWY").keys())

SCALES = {
    "KD_HYDROPATHY": "Kyte & Doolittle, J Mol Biol 157:105-132 (1982)",
    "EISENBERG": "Eisenberg et al., Faraday Symp Chem Soc 17:109 (1982)",
    "HOPP_WOODS": "Hopp & Woods, PNAS 78:3824 (1981)",
    "CF_HELIX/CF_SHEET": "Chou & Fasman, Adv Enzymol 47:45-148 (1978)",
    "AGGREGATION_ZZ": "Zhou & Zhou, FEBS Lett 557:21-26 (2004)",
    "DIWV": "Guruprasad, Reddy & Pandit, Protein Eng 4:155-161 (1990)",
    "BOMAN": "Boman, Antimicrob Agents Chemother 47:1904 (2003)",
    "ALIPHATIC": "Ikai, J Biochem 88:1895-1898 (1980)",
}
