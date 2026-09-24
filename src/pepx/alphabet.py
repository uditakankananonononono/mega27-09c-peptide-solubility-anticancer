"""Canonical amino-acid alphabet and validated peptide sequence container."""
from __future__ import annotations

from dataclasses import dataclass

AA_STANDARD = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_IDX = {aa: i for i, aa in enumerate(AA_STANDARD)}
IDX_TO_AA = {i: aa for aa, i in AA_TO_IDX.items()}


class InvalidSequenceError(ValueError):
    """Raised when a sequence contains residues outside the canonical 20."""


def validate_sequence(seq: str, allow_x: bool = False) -> str:
    """Return the uppercased sequence if valid, else raise InvalidSequenceError.

    Parameters
    ----------
    seq:
        Amino-acid sequence, any case, whitespace tolerated at the ends.
    allow_x:
        If True, the ambiguous residue 'X' is accepted and mapped to index 20.
    """
    s = seq.strip().upper()
    alphabet = set(AA_STANDARD) | ({"X"} if allow_x else set())
    bad = sorted(set(s) - alphabet)
    if bad:
        raise InvalidSequenceError(f"invalid residues {bad} in sequence of length {len(s)}")
    if not s:
        raise InvalidSequenceError("empty sequence")
    return s


def encode_indices(seq: str, allow_x: bool = False) -> list[int]:
    """Map a validated sequence to integer indices (X -> len(AA_STANDARD))."""
    s = validate_sequence(seq, allow_x=allow_x)
    x_idx = len(AA_STANDARD)
    return [AA_TO_IDX.get(aa, x_idx) for aa in s]


@dataclass(frozen=True)
class Peptide:
    """An immutable, validated peptide record."""

    sequence: str
    name: str = ""
    source: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "sequence", validate_sequence(self.sequence))

    def __len__(self) -> int:  # number of residues
        return len(self.sequence)

    @property
    def indices(self) -> list[int]:
        return encode_indices(self.sequence)
