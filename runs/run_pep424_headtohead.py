"""pep424 head-to-head: FoldAmyloid published predictions (AmyloGramAnalysis
benchmark files) vs my GNN/CNNv2 amyloid models scoring by hexapeptide
windowing. Labels: pep424_evaluation.txt (WALTZ-DB-derived)."""
import sys, json, time
sys.path.insert(0, 'src')
import numpy as np
import torch
from sklearn.metrics import roc_auc_score, matthews_corrcoef, accuracy_score
from pepx.datasets import load_amylogram
from pepx.models import PeptideGNN, PeptideCNNv2
from pepx.trainer import train_single_task, to_tensors, standardize

t0 = time.time()
RAW = 'data/raw/aggregation'
# --- labels ---
names, seqs, labels = [], [], []
for line in open(f'{RAW}/pep424_evaluation.txt'):
    parts = line.split()
    if len(parts) >= 3 and parts[-1] in ('+', '-'):
        names.append(parts[0]); seqs.append(parts[1]); labels.append(1 if parts[-1] == '+' else 0)
y = np.array(labels)
print(f"pep424 labels: {len(y)} ({y.sum()} pos)", flush=True)

# --- FoldAmyloid scores from published prediction file ---
fa_scores = {}
cur_name = None
for line in open(f'{RAW}/foldamyloid_pred.txt'):
    line = line.strip()
    if line.startswith('>'):
        cur_name = line[1:]
    elif line and set(line) <= {'0', '1'} and cur_name:
        bits = np.array([int(c) for c in line])
        fa_scores[cur_name] = bits.mean()  # fraction of residues in amyloid stretches
fa = np.array([fa_scores.get(n, np.nan) for n in names])
ok = ~np.isnan(fa)
fa_auc = roc_auc_score(y[ok], fa[ok])
print(f"FoldAmyloid on pep424: n={ok.sum()} AUC={fa_auc:.4f}", flush=True)

# --- my models: train on hexapeptides, score pep424 by window max-pool ---
full = load_amylogram('full')
bench = load_amylogram('benchmark')
bseqs = {r.sequence for r in bench}
tr = [r for r in full if r.sequence not in bseqs]

class Rec:
    def __init__(s, seq, label=0, meta=''): s.sequence, s.label, s.meta = seq, label, meta

from pepx.alphabet import AA_STANDARD as _AA_SET
def window_score(model, seq, max_len, wide, mu=None, sd=None):
    canon = set(_AA_SET)
    if not set(seq.upper()) <= canon:
        return np.nan  # non-canonical residues (e.g. O, U) - excluded, counted
    wins = [seq[i:i+6] for i in range(0, max(len(seq) - 5, 1))]
    wins = [w for w in wins if len(w) >= 6 and set(w) <= canon]
    if not wins:
        return np.nan
    recs = [Rec(w) for w in wins]
    i, d, _ = to_tensors(recs, 6, wide)
    if mu is not None:
        d = (d - mu) / sd
    with torch.no_grad():
        p = torch.sigmoid(model(i, d)).numpy()
    return float(p.max())

# GNN arm
res_gnn, _ = train_single_task(PeptideGNN(), tr, bench, 'amylogram', 'GNN-agg', max_len=32, epochs=20)
model = res_gnn_model = PeptideGNN()  # placeholder, retrained below via state
# NOTE: train_single_task returns metrics only; retrain capturing the model object
model = PeptideGNN()
res_gnn, _ = train_single_task(model, tr, bench, 'amylogram', 'GNN-agg', max_len=32, epochs=20)
print(res_gnn, flush=True)
my_scores = np.array([window_score(model, s, 6, False) for s in seqs])
ok2 = ~np.isnan(my_scores)
my_auc = roc_auc_score(y[ok2], my_scores[ok2])
print(f"pepx-GNN(hex-window) on pep424: n={ok2.sum()} AUC={my_auc:.4f} (excluded {(~ok2).sum()} non-canonical)", flush=True)

out = dict(dataset='pep424', n=int(len(y)), positives=int(y.sum()),
           foldamyloid_auc=float(fa_auc), pepx_gnn_auc=float(my_auc), n_scored=int(ok2.sum()), n_excluded_noncanonical=int((~ok2).sum()),
           beat_foldamyloid=bool(my_auc > fa_auc),
           source='FoldAmyloid predictions: michbur/AmyloGramAnalysis benchmark/FoldAmyloid_pred.txt; labels: pep424_evaluation.txt')
json.dump(out, open('results/pep424_headtohead.json', 'w'), indent=2)
print(out, flush=True)
print('elapsed', round(time.time() - t0, 1), 's', flush=True)
