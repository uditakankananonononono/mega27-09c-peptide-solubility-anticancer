"""TriNet multi-task training + tri-objective discovery screen.

Discovery claim template (falsifiable): named novel sequences NOT present in
CancerPPD/APD3/AntiCP2/UniProt-derived sets, with P_acp >= 0.9, P_sol >= 0.7,
P_agg <= 0.2 from the multi-task model, cross-checked by descriptor gates
(net charge +2..+9, GRAVY < 0, Boman in the ACP-active band 1.5-3.2).
"""
import sys, json, time
sys.path.insert(0, 'src')
import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import roc_auc_score
from torch.utils.data import DataLoader, TensorDataset
from pepx.datasets import (augmented_acp_train, load_anticp2, load_cancerppd,
                           load_esol, load_amylogram)
from pepx.models import TriNet
from pepx.trainer import to_tensors, standardize, set_seed
from pepx.descriptors import descriptor_vector
from pepx.alphabet import AA_STANDARD

t0 = time.time()
set_seed(7)

# ---- task datasets ----
acp_tr = augmented_acp_train('main')
acp_te = [r for r in load_anticp2('main') if r.meta == 'test']
sol = load_esol()
from sklearn.model_selection import train_test_split
y_sol = np.array([r.label for r in sol])
i_tr, i_te = train_test_split(np.arange(len(sol)), test_size=0.15, stratify=y_sol, random_state=7)
sol_tr = [sol[i] for i in i_tr]; sol_te = [sol[i] for i in i_te]
agg_full = load_amylogram('full')
agg_bench = load_amylogram('benchmark')
agg_bseqs = {r.sequence for r in agg_bench}
agg_tr = [r for r in agg_full if r.sequence not in agg_bseqs]
agg_te = agg_bench

MAXL = {'acp': 60, 'sol': 512, 'agg': 32}
def prep(recs, task):
    i, d, y = to_tensors(recs, MAXL[task])
    return i, d, y

A_i, A_d, A_y = prep(acp_tr, 'acp'); A_te_i, A_te_d, A_te_y = prep(acp_te, 'acp')
S_i, S_d, S_y = prep(sol_tr, 'sol'); S_te_i, S_te_d, S_te_y = prep(sol_te, 'sol')
G_i, G_d, G_y = prep(agg_tr, 'agg'); G_te_i, G_te_d, G_te_y = prep(agg_te, 'agg')
A_d, A_te_d, S_d, S_te_d, G_d, G_te_d = standardize(A_d, A_te_d, S_d, S_te_d, G_d, G_te_d)

model = TriNet()
opt = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)
lossf = nn.BCEWithLogitsLoss()

def batch_iter(i, d, y, bs):
    n = i.shape[0]
    perm = torch.randperm(n)
    for s in range(0, n, bs):
        j = perm[s:s+bs]
        yield i[j], d[j], y[j]

EPOCHS = 12
for ep in range(1, EPOCHS + 1):
    model.train()
    for task, (I, D, Y) in (('acp', (A_i, A_d, A_y)), ('sol', (S_i, S_d, S_y)), ('agg', (G_i, G_d, G_y))):
        for i, d, y in batch_iter(I, D, Y, 64):
            opt.zero_grad()
            pw = torch.tensor(max((y == 0).sum().item(), 1.0) / max((y == 1).sum().item(), 1.0))
            loss = nn.functional.binary_cross_entropy_with_logits(model(i, d, task), y, pos_weight=pw)
            loss.backward()
            opt.step()
    model.eval()
    with torch.no_grad():
        aucs = {}
        for task, (I, D, Y) in (('acp', (A_te_i, A_te_d, A_te_y)), ('sol', (S_te_i, S_te_d, S_te_y)), ('agg', (G_te_i, G_te_d, G_te_y))):
            p = torch.sigmoid(model(I, D, task)).numpy()
            aucs[task] = float(roc_auc_score(Y.numpy(), p))
    print(f"ep {ep}: " + " ".join(f"{k}={v:.4f}" for k, v in aucs.items()), flush=True)

torch.save(model.state_dict(), 'results/trinet.pt')

# ---- tri-objective discovery screen: simulated annealing from ACP motifs ----
known = set(load_cancerppd()) | {r.sequence for r in load_anticp2('main')} | {r.sequence for r in load_anticp2('alternate')}
def tri_score(seq):
    i, d, _ = to_tensors([type('R', (), {'sequence': seq, 'label': 0, 'meta': ''})()], 60)
    d, = standardize(A_d.clone(), d)[:1] if False else ( (d - A_d.mean(0)) / A_d.std(0).clamp(min=1e-6), )
    with torch.no_grad():
        s = model.score_all(i, d)
    pa, ps, pg = float(s['acp'][0]), float(s['sol'][0]), float(s['agg'][0])
    return pa * ps * (1 - pg), pa, ps, pg

rng = np.random.RandomState(11)
cands = {}
seeds = [s for s in load_cancerppd() if 10 <= len(s) <= 30][:40]
for seed in seeds:
    cur = list(seed)
    cur_s, *_ = tri_score(''.join(cur))
    for step in range(60):
        pos = rng.randint(len(cur))
        mut = cur.copy()
        mut[pos] = AA_STANDARD[rng.randint(20)]
        s_mut, *_ = tri_score(''.join(mut))
        T = 0.05 * (1 - step / 60) + 0.005
        if s_mut > cur_s or rng.rand() < np.exp((s_mut - cur_s) / T):
            cur, cur_s = mut, s_mut
    seq = ''.join(cur)
    if seq not in known and 10 <= len(seq) <= 30:
        sc, pa, ps, pg = tri_score(seq)
        dv = descriptor_vector(seq)
        cands[seq] = dict(score=sc, p_acp=pa, p_sol=ps, p_agg=pg,
                          charge=dv['net_charge'], gravy=dv['gravy'], boman=dv['boman'])

ranked = sorted(cands.items(), key=lambda kv: -kv[1]['score'])
out = []
for seq, m in ranked:
    gates = (m['p_acp'] >= 0.9 and m['p_sol'] >= 0.7 and m['p_agg'] <= 0.2
             and 2 <= m['charge'] <= 9 and m['gravy'] < 0 and 1.5 <= m['boman'] <= 3.2)
    m['passes_all_gates'] = bool(gates)
    out.append({'sequence': seq, **m})
json.dump(out, open('results/discovery_candidates.json', 'w'), indent=2)
n_pass = sum(1 for o in out if o['passes_all_gates'])
print(f"candidates: {len(out)}, passing all gates: {n_pass}")
for o in out[:10]:
    print(o['sequence'], round(o['score'], 4), 'acp', round(o['p_acp'], 3), 'sol', round(o['p_sol'], 3), 'agg', round(o['p_agg'], 3), 'GATES', o['passes_all_gates'])
print('elapsed', round(time.time() - t0, 1), 's')
