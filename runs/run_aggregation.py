"""Amyloid aggregation benchmark: train on AmyloGram full, locked test on the
AmyloGram benchmark split, decontaminated by exact sequence."""
import sys, json, time
sys.path.insert(0, 'src')
import numpy as np
from pepx.datasets import load_amylogram
from pepx.baselines import run_benchmark
from pepx.models import PeptideCNNv2, PeptideGNN
from pepx.trainer import train_single_task

t0 = time.time()
full = load_amylogram('full')
bench = load_amylogram('benchmark')
bench_seqs = {r.sequence for r in bench}
train = [r for r in full if r.sequence not in bench_seqs]
test = bench
print(f"train={len(train)} test={len(test)} (decontaminated)")

for r in run_benchmark(train, test, 'amylogram'):
    print(f"{r.model:8s} acc={r.accuracy:.4f} mcc={r.mcc:.4f} auc={r.auc:.4f}")

res, _ = train_single_task(PeptideCNNv2(), train, test, 'amylogram', 'CNNv2',
                           max_len=32, epochs=30, wide=True)
print(res)
res2, _ = train_single_task(PeptideGNN(), train, test, 'amylogram', 'GNN',
                            max_len=32, epochs=30)
print(res2)
print('elapsed', round(time.time() - t0, 1), 's')
