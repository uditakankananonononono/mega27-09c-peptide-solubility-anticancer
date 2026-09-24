"""Minimal, strict FASTA reader/writer used by all dataset loaders."""
from __future__ import annotations

from pathlib import Path
from typing import Iterator, Tuple


def read_fasta(path: str | Path) -> Iterator[Tuple[str, str]]:
    """Yield (header, sequence) pairs. Sequences are concatenated across lines.

    Raises ValueError on malformed input (sequence before first header).
    """
    header: str | None = None
    chunks: list[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                if header is None:
                    raise ValueError(f"{path}:{lineno}: sequence data before first header")
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def write_fasta(records: Iterator[Tuple[str, str]], path: str | Path, wrap: int = 80) -> int:
    """Write (header, sequence) pairs to FASTA; returns record count."""
    n = 0
    with open(path, "w", encoding="utf-8") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), wrap):
                fh.write(seq[i : i + wrap] + "\n")
            n += 1
    return n
