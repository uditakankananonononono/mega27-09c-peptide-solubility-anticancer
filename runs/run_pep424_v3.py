"""pep424 head-to-head v3 - ORDER-ALIGNED join (names are not unique:
158 entries share the header '3D profile'). pep424.fasta (419) and
FoldAmyloid_pred.txt (419) are in the same order; labels from
pep424_evaluation.txt joined by SEQUENCE.
(a) FoldAmyloid published-prediction AUC vs (b) pepx-GNN transfer AUC on the
identical labeled subset; (c) 5-fold CV vs AmyloGram published CV 0.865."""
import sys, json, time
sys.path.insert(0, 'src')
import numpy as np
import torch
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from pepx.datasets import load_amylogram
from pepx.fasta import read_fasta
from pepx.models import PeptideGNN
from pepx.trainer import train_single_task, to_tensors
from pepx.encoders import descriptor_array

t0 = time.time()
RAW = 'data/raw/aggregation'
seq2label = {}
for line in open(f'{RAW}/pep424_evaluation.txt'):
    parts = line.rstrip('\n').split('\t')
    if len(parts) >= 3 and parts[2].strip() in ('+', '-'):
        seq2label[parts[1].strip().upper()] = 1 if parts[2].strip() == '+' else 0

fasta_seqs = [s.upper() for _, s in read_fasta(f'{RAW}/pep424.fasta')]
assert len(fasta_seqs) == 419, len(fasta_seqs)

fa_scores = []
cur_bits = None
for line in open(f'{RAW}/foldamyloid_pred.txt'):
    line = line.strip()
    if line.startswith('>'):
        if cur_bits is not None:
            fa_scores.append(cur_bits)
        cur_bits = ''
    elif line and set(line) <= {'0', '1'}:
        cur_bits = (cur_bits or '') + line
if cur_bits:
    fa_scores.append(cur_bits)
fa = np.array([np.array([int(c) for c in b]).mean() if b else np.nan for b in fa_scores])
assert len(fa) == 419, len(fa)

y_all = np.array([seq2label.get(s, -1) for s in fasta_seqs])
ok = y_all >= 0
y = y_all[ok]
fa_ok = fa[ok]
seqs = [s for s, k in zip(fasta_seqs, ok) if k]
print(f"aligned: {ok.sum()} of 419 with labels ({y.sum()} pos)", flush=True)
fa_auc = roc_auc_score(y, fa_ok)
print(f"FoldAmyloid published preds: AUC={fa_auc:.4f}", flush=True)

class Rec:
    def __init__(s, seq, label=0, meta=''):
        s.sequence, s.label, s.meta = seq, label, meta

full = load_amylogram('full')
bench = load_amylogram('benchmark')
bseqs = {r.sequence for r in bench}
tr = [r for r in full if r.sequence not in bseqs]
Dtr = np.stack([descriptor_array(r.sequence) for r in tr])
mu, sd = Dtr.mean(0), Dtr.std(0); sd[sd < 1e-6] = 1e-6
model = PeptideGNN()
res_gnn, _ = train_single_task(model, tr, bench, 'amylogram', 'GNN-agg', max_len=32, epochs=20)
print(res_gnn, flush=True)

def score_seqs(model, sequences):
    model.eval(); out = []
    for s in range(0, len(sequences), 128):
        recs = [Rec(x) for x in sequences[s:s+128]]
        i, d, _ = to_tensors(recs, 32)
        d = (d - torch.from_numpy(mu).float()) / torch.from_numpy(sd).float()
        with torch.no_grad():
            out.append(torch.sigmoid(model(i, d)).numpy())
    return np.concatenate(out)

p_tr = score_seqs(model, seqs)
transfer_auc = roc_auc_score(y, p_tr)
print(f"TRANSFER pepx-GNN AUC={transfer_auc:.4f} vs FoldAmyloid {fa_auc:.4f}", flush=True)

recs = [Rec(s, int(l)) for s, l in zip(seqs, y)]
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=7)
cv_prob = np.zeros(len(recs))
for fold, (itr, ite) in enumerate(skf.split(np.zeros(len(recs)), y)):
    m = PeptideGNN()
    _, p = train_single_task(m, [recs[i] for i in itr], [recs[i] for i in ite],
                             f'pep424-cv{fold}', 'GNN', max_len=32, epochs=20, seed=7 + fold)
    cv_prob[ite] = p
cv_auc = roc_auc_score(y, cv_prob)
print(f"CV5 pepx-GNN AUC={cv_auc:.4f} vs AmyloGram published 0.865", flush=True)

out = dict(dataset='pep424', aligned_n=int(ok.sum()), positives=int(y.sum()),
           foldamyloid_auc=float(fa_auc), transfer_pepx_gnn_auc=float(transfer_auc),
           transfer_beats_foldamyloid=bool(transfer_auc > fa_auc),
           cv5_pepx_gnn_auc=float(cv_auc), amylogram_published_cv_auc=0.865,
           join='order-aligned fasta<->foldamyloid preds; labels by sequence',
           sources='michbur/AmyloGramAnalysis benchmark/{pep424.fasta,FoldAmyloid_pred.txt,pep424_evaluation.txt}')
json.dump(out, open('results/pep424_v3.json', 'w'), indent=2)
print(out, flush=True)
print('elapsed', round(time.time() - t0, 1), 's', flush=True)
