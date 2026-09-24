"""Dataset loaders. Every loader reads from data/raw and returns clean,
validated (sequence, label, meta) records. Provenance is logged to
data/MANIFEST.md - no silent data."""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

from .alphabet import AA_STANDARD, InvalidSequenceError, validate_sequence
from .fasta import read_fasta

RAW = Path(__file__).resolve().parents[2] / "data" / "raw"


@dataclass(frozen=True)
class LabeledPeptide:
    sequence: str
    label: int  # 1 = positive class, 0 = negative
    source: str
    meta: str = ""


def _clean_lines(path: Path) -> List[str]:
    out = []
    for line in path.read_text().splitlines():
        s = line.strip().upper()
        if s and set(s) <= set(AA_STANDARD):
            out.append(s)
    return out


def load_anticp2(split: str = "main") -> List[LabeledPeptide]:
    """AntiCP 2.0 benchmark (Agrawal et al. 2021): main = 689+689 train /
    172+172 test experimentally validated ACPs vs non-ACP AMPs; alternate =
    777+777 / 193+193 vs random UniProt-derived negatives."""
    base = RAW / "anticp2"
    recs: List[LabeledPeptide] = []
    for part in ("train", "test"):
        for cls, lab in (("pos", 1), ("neg", 0)):
            f = base / f"{cls}_{part}_{split}"
            for s in _clean_lines(f):
                recs.append(LabeledPeptide(s, lab, f"anticp2_{split}", part))
    return recs


def load_cancerppd(min_len: int = 5, max_len: int = 100) -> List[str]:
    """CancerPPD (Tyagi et al. 2015): experimentally validated anticancer
    peptides, L-natural canonical subset."""
    seqs: List[str] = []
    with open(RAW / "cancerppd" / "l_natural.txt") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 2:
                continue
            try:
                s = validate_sequence(row[1])
            except InvalidSequenceError:
                continue
            if min_len <= len(s) <= max_len:
                seqs.append(s)
    return seqs


def load_esol(threshold_pct: float = 30.0) -> List[LabeledPeptide]:
    """eSOL E. coli K-12 solubility (Niwa et al. PNAS 2009; Nucleic Acids Res
    2012): PURE-system measured solubility % for 4,132 proteins, sequences
    mapped from UniProt proteome UP000000625 by gene name. Binary label at the
    published 30% cutoff used by solubility-classification benchmarks."""
    seqs: Dict[str, str] = {}
    for header, s in read_fasta(RAW / "esol" / "ecoli_k12_proteome.fasta"):
        gene = None
        for tok in header.split():
            if tok.startswith("GN="):
                gene = tok[3:]
                break
        if gene and gene not in seqs:
            seqs[gene] = s
    recs: List[LabeledPeptide] = []
    with open(RAW / "esol" / "esol.csv") as fh:
        for row in csv.DictReader(fh):
            gene = (row.get("Gene name K-12") or "").strip()
            sol = (row.get("Solubility (%)") or "").strip()
            if not gene or not sol or gene not in seqs:
                continue
            try:
                pct = float(sol)
                s = validate_sequence(seqs[gene])
            except (ValueError, InvalidSequenceError):
                continue
            recs.append(LabeledPeptide(s, int(pct >= threshold_pct), "esol", f"solubility_pct={pct}"))
    return recs


def load_amylogram(part: str = "full") -> List[LabeledPeptide]:
    """AmyloGram/WALTZ-DB-derived amyloid hexapeptide benchmark (Kozlowski &
    Burdukiewicz 2017; WALTZ-DB 2.0, Beerten et al. 2015). part: 'full'
    (421+1044) or 'benchmark' (269+746)."""
    base = RAW / "aggregation"
    suffix = "full" if part == "full" else "benchmark"
    recs: List[LabeledPeptide] = []
    for cls, lab in (("pos", 1), ("neg", 0)):
        for header, s in read_fasta(base / f"amyloid_{cls}_{suffix}.fasta"):
            try:
                recs.append(LabeledPeptide(validate_sequence(s), lab, f"amylogram_{suffix}", header))
            except InvalidSequenceError:
                continue
    return recs


def load_apd3(max_len: int = 100) -> List[str]:
    """APD3 natural AMPs (Wang et al. 2024 update) - non-anticancer-negative
    augmentation pool after decontamination against CancerPPD and AntiCP2."""
    seqs = []
    for header, s in read_fasta(RAW / "apd3" / "naturalAMPs_APD2024a.fasta"):
        try:
            v = validate_sequence(s)
            if 5 <= len(v) <= max_len:
                seqs.append(v)
        except InvalidSequenceError:
            continue
    return seqs


def augmented_acp_train(split: str = "main") -> List[LabeledPeptide]:
    """AntiCP2 train split + CancerPPD ACPs (pos) + APD3 AMPs (neg), fully
    decontaminated against the locked AntiCP2 test split and each other."""
    recs = load_anticp2(split)
    train = [r for r in recs if r.meta == "train"]
    test_seqs = {r.sequence for r in recs if r.meta == "test"}
    train_seqs = {r.sequence for r in train}
    acps = set(load_cancerppd())
    out = list(train)
    for s in sorted(acps - test_seqs - train_seqs):
        out.append(LabeledPeptide(s, 1, "cancerppd_aug"))
    pos_all = acps | {r.sequence for r in train if r.label == 1}
    for s in sorted(set(load_apd3()) - pos_all - test_seqs - train_seqs):
        out.append(LabeledPeptide(s, 0, "apd3_aug"))
    return out
