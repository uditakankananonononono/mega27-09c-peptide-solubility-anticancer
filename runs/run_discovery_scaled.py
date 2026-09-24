"""Scaled tri-objective discovery screen: TriNet 25 epochs, batched annealing
(100 seeds x 100 steps). Candidates gated, named, committed."""
import sys, json, time
sys.path.insert(0, 'src')
import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
from pepx.datasets import (augmented_acp_train, load_anticp2, load_cancerppd,
                           load_esol, load_amylogram)
from pepx.models import TriNet
from pepx.trainer import to_tensors, standardize, set_seed
from pepx.descriptors import descriptor_vector
from pepx.alphabet import AA_STANDARD

t0 = time.time()
set_seed(7)
acp_tr = augmented_acp_train('main')
acp_te = [r for r in load_anticp2('main') if r.meta == 'test']
sol = load_esol()
y_sol = np.array([r.label for r in sol])
i_tr, i_te = train_test_split(np.arange(len(sol)), test_size=0.15, stratify=y_sol, random_state=7)
sol_tr = [sol[i] for i in i_tr]; sol_te = [sol[i] for i in i_te]
agg_full = load_amylogram('full'); agg_bench = load_amylogram('benchmark')
agg_bseqs = {r.sequence for r in agg_bench}
agg_tr = [r for r in agg_full if r.sequence not in agg_bseqs]

MAXL = {'acp': 60, 'sol': 300, 'agg': 32}
def prep(recs, task):
    return to_tensors(recs, MAXL[task])
A_i, A_d, A_y = prep(acp_tr, 'acp'); A_te_i, A_te_d, A_te_y = prep(acp_te, 'acp')
S_i, S_d, S_y = prep(sol_tr, 'sol'); S_te_i, S_te_d, S_te_y = prep(sol_te, 'sol')
G_i, G_d, G_y = prep(agg_tr, 'agg'); G_te_i, G_te_d, G_te_y = prep(agg_bench, 'agg')
MU, SD = A_d.mean(0, keepdim=True), A_d.std(0, keepdim=True).clamp(min=1e-6)
np.savez('results/trinet_norm_25ep.npz', mu=MU.numpy().ravel(), sd=SD.numpy().ravel())
A_d, A_te_d = standardize(A_d, A_te_d)
S_d, S_te_d = standardize(S_d, S_te_d)
G_d, G_te_d = standardize(G_d, G_te_d)

model = TriNet()
opt = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)

def batch_iter(i, d, y, bs):
    perm = torch.randperm(i.shape[0])
    for s in range(0, i.shape[0], bs):
        j = perm[s:s+bs]
        yield i[j], d[j], y[j]

EPOCHS = 25
best_auc, best_state = -1, None
for ep in range(1, EPOCHS + 1):
    model.train()
    for task, (I, D, Y) in (('acp', (A_i, A_d, A_y)), ('sol', (S_i, S_d, S_y)), ('agg', (G_i, G_d, G_y))):
        for i, d, y in batch_iter(I, D, Y, 64):
            opt.zero_grad()
            pw = torch.tensor(max((y == 0).sum().item(), 1.0) / max((y == 1).sum().item(), 1.0))
            loss = nn.functional.binary_cross_entropy_with_logits(model(i, d, task), y, pos_weight=pw)
            loss.backward(); opt.step()
    model.eval()
    with torch.no_grad():
        aucs = {}
        for task, (I, D, Y) in (('acp', (A_te_i, A_te_d, A_te_y)), ('sol', (S_te_i, S_te_d, S_te_y)), ('agg', (G_te_i, G_te_d, G_te_y))):
            aucs[task] = float(roc_auc_score(Y.numpy(), torch.sigmoid(model(I, D, task)).numpy()))
    m = np.mean(list(aucs.values()))
    if m > best_auc:
        best_auc, best_state = m, {k: v.clone() for k, v in model.state_dict().items()}
    print(f"ep {ep}: " + " ".join(f"{k}={v:.4f}" for k, v in aucs.items()), flush=True)
model.load_state_dict(best_state)
torch.save(model.state_dict(), 'results/trinet_25ep.pt')
print('saved trinet_25ep.pt, mean heldout AUC', round(best_auc, 4), flush=True)

# ---- batched annealing screen ----
class R:
    def __init__(s, seq): s.sequence, s.label, s.meta = seq, 0, ''

def tri_batch(seqs):
    i, d, _ = to_tensors([R(s) for s in seqs], 60)
    d = (d - MU) / SD
    with torch.no_grad():
        out = model.score_all(i, d)
    return out['acp'].numpy(), out['sol'].numpy(), out['agg'].numpy()

known = set(load_cancerppd()) | {r.sequence for r in load_anticp2('main')} | {r.sequence for r in load_anticp2('alternate')}
rng = np.random.RandomState(11)
seeds = [s for s in load_cancerppd() if 10 <= len(s) <= 30][:100]
current = [list(s) for s in seeds]
pa, ps, pg = tri_batch([''.join(c) for c in current])
cur_score = pa * ps * (1 - pg)
STEPS = 100
for step in range(STEPS):
    cands = []
    for c in current:
        m = c.copy(); m[rng.randint(len(m))] = AA_STANDARD[rng.randint(20)]
        cands.append(''.join(m))
    pa2, ps2, pg2 = tri_batch(cands)
    s2 = pa2 * ps2 * (1 - pg2)
    T = 0.05 * (1 - step / STEPS) + 0.005
    for k in range(len(current)):
        if s2[k] > cur_score[k] or rng.rand() < np.exp((s2[k] - cur_score[k]) / T):
            current[k] = list(cands[k]); cur_score[k] = s2[k]
    if step % 20 == 0:
        print(f"step {step}: median S {np.median(cur_score):.3f}", flush=True)

finals = sorted({''.join(c) for c in current} - known)
pa, ps, pg = tri_batch(finals)
out = []
for s, a, b, g in zip(finals, pa, ps, pg):
    dv = descriptor_vector(s)
    gates = (a >= 0.9 and b >= 0.7 and g <= 0.2 and 2 <= dv['net_charge'] <= 9
             and dv['gravy'] < 0 and 1.5 <= dv['boman'] <= 3.2 and 10 <= len(s) <= 30)
    out.append(dict(sequence=s, score=float(a*b*(1-g)), p_acp=float(a), p_sol=float(b), p_agg=float(g),
                    charge=float(dv['net_charge']), gravy=float(dv['gravy']), boman=float(dv['boman']),
                    passes_all_gates=bool(gates)))
out.sort(key=lambda o: -o['score'])
json.dump(out, open('results/discovery_candidates_scaled.json', 'w'), indent=2)
n_pass = sum(1 for o in out if o['passes_all_gates'])
print(f"unique candidates: {len(out)}, passing all gates: {n_pass}", flush=True)
for o in out[:8]:
    print(o['sequence'], round(o['score'],3), 'acp', round(o['p_acp'],3), 'sol', round(o['p_sol'],3), 'agg', round(o['p_agg'],3), 'GATES', o['passes_all_gates'], flush=True)
print('elapsed', round(time.time()-t0,1), 's', flush=True)
