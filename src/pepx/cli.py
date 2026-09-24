"""pepx - tri-objective peptide intelligence CLI.

Usage:
  python -m pepx.cli score SEQUENCE [SEQUENCE ...]
  python -m pepx.cli score-file INPUT.fasta [--out results.tsv]
  python -m pepx.cli describe SEQUENCE
Scores every peptide on anticancer probability, solubility probability and
aggregation probability (TriNet multi-task model) plus the full published
descriptor panel. Ships with the trained weights in results/trinet.pt.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import torch

from .alphabet import InvalidSequenceError, validate_sequence
from .descriptors import descriptor_vector
from .fasta import read_fasta
from .models import TriNet
from .trainer import to_tensors

ROOT = Path(__file__).resolve().parents[2]


class _Rec:
    def __init__(self, seq):
        self.sequence, self.label, self.meta = seq, 0, ""


def _load_trinet(path: Path):
    m = TriNet()
    m.load_state_dict(torch.load(path, map_location="cpu"))
    m.eval()
    norm = path.with_name("trinet_norm.npz")
    mu = sd = None
    if norm.exists():
        z = np.load(norm)
        mu, sd = z["mu"], z["sd"]
    return m, mu, sd


def score_sequences(model: TriNet, seqs, mu=None, sd=None):
    recs = [_Rec(s) for s in seqs]
    i, d, _ = to_tensors(recs, 60)
    if mu is not None:
        d = (d - torch.from_numpy(mu).float()) / torch.from_numpy(sd).float()
    with torch.no_grad():
        s = model.score_all(i, d)
    return {t: v.numpy() for t, v in s.items()}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="pepx", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("score", help="score sequences on acp/sol/agg")
    p1.add_argument("sequences", nargs="+")
    p1.add_argument("--model", default=str(ROOT / "results" / "trinet.pt"))
    p2 = sub.add_parser("score-file", help="score a FASTA file")
    p2.add_argument("fasta")
    p2.add_argument("--out", default=None)
    p2.add_argument("--model", default=str(ROOT / "results" / "trinet.pt"))
    p3 = sub.add_parser("describe", help="physicochemical descriptor panel")
    p3.add_argument("sequence")
    args = ap.parse_args(argv)

    if args.cmd == "describe":
        v = descriptor_vector(validate_sequence(args.sequence))
        print(json.dumps(v, indent=2))
        return 0

    model_path = Path(args.model)
    if not model_path.exists():
        print(f"error: trained weights not found at {model_path}", file=sys.stderr)
        return 2
    model, mu, sd = _load_trinet(model_path)

    if args.cmd == "score":
        seqs = [validate_sequence(s) for s in args.sequences]
        scores = score_sequences(model, seqs, mu, sd)
        print("sequence\tP_anticancer\tP_soluble\tP_aggregation")
        for j, s in enumerate(seqs):
            print(f"{s}\t{scores['acp'][j]:.4f}\t{scores['sol'][j]:.4f}\t{scores['agg'][j]:.4f}")
        return 0

    if args.cmd == "score-file":
        rows = list(read_fasta(args.fasta))
        seqs, keep = [], []
        for h, s in rows:
            try:
                seqs.append(validate_sequence(s)); keep.append(h)
            except InvalidSequenceError:
                print(f"skip invalid: {h}", file=sys.stderr)
        scores = score_sequences(model, seqs, mu, sd)
        lines = ["name\tsequence\tP_anticancer\tP_soluble\tP_aggregation"]
        for j, h in enumerate(keep):
            lines.append(f"{h}\t{seqs[j]}\t{scores['acp'][j]:.4f}\t{scores['sol'][j]:.4f}\t{scores['agg'][j]:.4f}")
        out = "\n".join(lines)
        if args.out:
            Path(args.out).write_text(out + "\n")
            print(f"wrote {len(keep)} scores to {args.out}")
        else:
            print(out)
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
