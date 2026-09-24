"""Sequence -> numeric feature encoders shared by all models."""
from __future__ import annotations

import numpy as np

from .alphabet import AA_STANDARD, AA_TO_IDX, validate_sequence
from .descriptors import DESCRIPTOR_NAMES, descriptor_vector

N_AA = len(AA_STANDARD)


def one_hot(seq: str, max_len: int) -> np.ndarray:
    """(max_len, 21) one-hot; index 20 = X/pad-mask channel is zeros for pad."""
    s = validate_sequence(seq)
    out = np.zeros((max_len, N_AA + 1), dtype=np.float32)
    for i, aa in enumerate(s[:max_len]):
        out[i, AA_TO_IDX.get(aa, N_AA)] = 1.0
    return out


def descriptor_array(seq: str) -> np.ndarray:
    v = descriptor_vector(seq)
    return np.array([v[k] for k in DESCRIPTOR_NAMES], dtype=np.float32)


def kmer_counts(seq: str, k: int = 2) -> np.ndarray:
    """Normalized k-mer composition vector of length 20**k (k<=3)."""
    s = validate_sequence(seq)
    dim = N_AA ** k
    out = np.zeros(dim, dtype=np.float32)
    idx = [AA_TO_IDX[a] for a in s]
    n = 0
    for i in range(len(idx) - k + 1):
        j = 0
        for off in range(k):
            j = j * N_AA + idx[i + off]
        out[j] += 1.0
        n += 1
    return out / max(n, 1)


def aac(seq: str) -> np.ndarray:
    """Amino-acid composition (20,)."""
    s = validate_sequence(seq)
    out = np.zeros(N_AA, dtype=np.float32)
    for a in s:
        out[AA_TO_IDX[a]] += 1.0
    return out / len(s)


def dpc(seq: str) -> np.ndarray:
    """Dipeptide composition (400,) - the AntiCP baseline encoding."""
    return kmer_counts(seq, 2)
