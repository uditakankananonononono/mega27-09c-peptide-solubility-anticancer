"""pep424 head-to-head v2 (fixed parse + fixed standardization):
(a) TRANSFER: GNN trained on AmyloGram hexapeptides scores pep424 directly,
    descriptors standardized with TRAIN stats. Comparison: FoldAmyloid's
    published predictions on the identical set (no pep424 training either).
(b) CV: 5-fold CV on pep424 (AmyloGram paper protocol) vs their published
    cross-validated AUC 0.865 (Kozlowski & Burdukiewicz 2017)."""
import sys, json, time
sys.path.insert(0, 'src')
import numpy as np
import torch
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from pepx.datasets import load_amylogram, LabeledPeptide
from pepx.models import PeptideGNN
from pepx.trainer import train_single_task, to_tensors, set_seed
from pepx.encoders import descriptor_array

t0 = time.time()
RAW = 'data/raw/aggregation'
names, seqs, labels = [], [], []
for line in open(f'{RAW}/pep424_evaluation.txt'):
    parts = line.rstrip('\n').split('\t')
    if len(parts) >= 3 and parts[2].strip() in ('+', '-'):
        names.append(parts[0]); seqs.append(parts[1].strip().upper()); labels.append(1 if parts[2].strip() == '+' else 0)
y = np.array(labels)
print(f"pep424: {len(y)} sequences, {y.sum()} positives", flush=True)

fa_scores = {}
cur = None
for line in open(f'{RAW}/foldamyloid_pred.txt'):
    line = line.strip()
    if line.startswith('>'):
        cur = line[1:]
    elif line and set(line) <= {'0', '1'} and cur:
        fa_scores[cur] = np.array([int(c) for c in line]).mean()
fa = np.array([fa_scores.get(n, np.nan) for n in names])
okfa = ~np.isnan(fa)
fa_auc = roc_auc_score(y[okfa], fa[okfa])
print(f"FoldAmyloid (published preds): n={okfa.sum()} AUC={fa_auc:.4f}", flush=True)

class Rec:
    def __init__(s, seq, label=0, meta=''):
        s.sequence, s.label, s.meta = seq, label, meta

# ---- (a) transfer ----
full = load_amylogram('full')
bench = load_amylogram('benchmark')
bseqs = {r.sequence for r in bench}
tr = [r for r in full if r.sequence not in bseqs]
mu = np.stack([descriptor_array(r.sequence) for r in tr]).mean(0)
sd = np.stack([descriptor_array(r.sequence) for r in tr]).std(0); sd[sd < 1e-6] = 1e-6

model = PeptideGNN()
res_gnn, _ = train_single_task(model, tr, bench, 'amylogram', 'GNN-agg', max_len=32, epochs=20)
print(res_gnn, flush=True)

def score_seqs(model, sequences, mu, sd):
    model.eval()
    out = []
    B = 128
    for s in range(0, len(sequences), B):
        recs = [Rec(x) for x in sequences[s:s+B]]
        i, d, _ = to_tensors(recs, 32)
        d = (d - torch.from_numpy(mu).float()) / torch.from_numpy(sd).float()
        with torch.no_grad():
            out.append(torch.sigmoid(model(i, d)).numpy())
    return np.concatenate(out)

p_transfer = score_seqs(model, seqs, mu, sd)
transfer_auc = roc_auc_score(y, p_transfer)
print(f"TRANSFER pepx-GNN: AUC={transfer_auc:.4f} (bar: FoldAmyloid {fa_auc:.4f})", flush=True)

# ---- (b) 5-fold CV on pep424 ----
recs = [Rec(s, int(l)) for s, l in zip(seqs, labels)]
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=7)
cv_prob = np.zeros(len(recs))
for fold, (itr, ite) in enumerate(skf.split(np.zeros(len(recs)), y)):
    tr_r = [recs[i] for i in itr]; te_r = [recs[i] for i in ite]
    m = PeptideGNN()
    _, p = train_single_task(m, tr_r, te_r, f'pep424-cv{fold}', 'GNN', max_len=32, epochs=20, seed=7 + fold)
    cv_prob[ite] = p
    print(f"fold {fold} done", flush=True)
cv_auc = roc_auc_score(y, cv_prob)
print(f"CV pepx-GNN on pep424: AUC={cv_auc:.4f} (AmyloGram published CV ~0.865)", flush=True)

out = dict(dataset='pep424', n=int(len(y)), positives=int(y.sum()),
           foldamyloid_auc=float(fa_auc),
           transfer_pepx_gnn_auc=float(transfer_auc),
           transfer_beats_foldamyloid=bool(transfer_auc > fa_auc),
           cv5_pepx_gnn_auc=float(cv_auc),
           amylogram_published_cv_auc=0.865,
           protocol='transfer: train on AmyloGram hexapeptides, score pep424; cv: 5-fold stratified seed 7')
json.dump(out, open('results/pep424_v2.json', 'w'), indent=2)
print(out, flush=True)
print('elapsed', round(time.time() - t0, 1), 's', flush=True)
