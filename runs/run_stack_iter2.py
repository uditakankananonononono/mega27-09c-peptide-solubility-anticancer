"""Iteration-2 per PREREG_ADDENDUM_2026-09-27 (locked): a priori equal-weight
ensembles; combination rule selected by train-internal 5-fold CV; single
official-test score; 10k bootstrap vs published 0.83/0.51."""
import sys, json, time
sys.path.insert(0, 'src')
t0 = time.time()
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, matthews_corrcoef
import torch
from pepx.datasets import load_anticp2
from pepx.encoders import aac, dpc
from pepx.models import PeptideCNNv2
from pepx.trainer import to_tensors, standardize, set_seed

SEED = 2709
recs = load_anticp2('main')
train = [r for r in recs if r.meta == 'train']
test = [r for r in recs if r.meta == 'test']
ytr = np.array([r.label for r in train]); yte = np.array([r.label for r in test])
MAXLEN = 50
Hp = dict(emb_dim=32, ch=96, blocks=2, dropout=0.2780616350232075)
LR = 0.0017286476459924243; EPOCHS = 18; NSEEDS = 7

def cnn_once(tr, ev, seed):
    set_seed(seed)
    m = PeptideCNNv2(**Hp)
    Xi, Xd, y = to_tensors(tr, MAXLEN, wide=True)
    Ei, Ed, _ = to_tensors(ev, MAXLEN, wide=True)
    Xd, Ed = standardize(Xd, Ed)
    opt = torch.optim.Adam(m.parameters(), lr=LR)
    lf = torch.nn.BCEWithLogitsLoss(); n = len(y)
    for ep in range(EPOCHS):
        m.train(); perm = torch.randperm(n)
        for i in range(0, n, 64):
            b = perm[i:i+64]; opt.zero_grad()
            lf(m(Xi[b], Xd[b]), y[b]).backward(); opt.step()
    m.eval()
    with torch.no_grad():
        return torch.sigmoid(m(Ei, Ed)).numpy()

def cnn_ens(tr, ev, base_seed):
    return np.mean([cnn_once(tr, ev, base_seed + s) for s in range(NSEEDS)], axis=0)

def dpc_arr(rs): return np.stack([dpc(r.sequence) for r in rs])
def aac_arr(rs): return np.stack([aac(r.sequence) for r in rs])

def fit_et(tr, y): 
    et = ExtraTreesClassifier(n_estimators=400, random_state=SEED, n_jobs=2)
    et.fit(dpc_arr(tr), y); return et
def fit_lg(tr, y):
    lg = LogisticRegression(max_iter=2000); lg.fit(aac_arr(tr), y); return lg

# --- selection CV: OOF probs for each base ---
oof_cnn = np.zeros(len(train)); oof_et = np.zeros(len(train)); oof_lg = np.zeros(len(train))
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
for fold, (ti, vi) in enumerate(skf.split(train, ytr)):
    tr_r = [train[i] for i in ti]; va_r = [train[i] for i in vi]; y_ti = ytr[ti]
    oof_cnn[vi] = cnn_ens(tr_r, va_r, SEED + 1000 * fold)
    oof_et[vi] = fit_et(tr_r, y_ti).predict_proba(dpc_arr(va_r))[:, 1]
    oof_lg[vi] = fit_lg(tr_r, y_ti).predict_proba(aac_arr(va_r))[:, 1]
    print('fold', fold, 'done', round(time.time() - t0), 's', flush=True)
    np.savez('/tmp/iter2_oof.npz', oof_cnn=oof_cnn, oof_et=oof_et, oof_lg=oof_lg, ytr=ytr)

cand = {
    'i_cnn7': oof_cnn,
    'ii_cnn7_et': 0.5 * oof_cnn + 0.5 * oof_et,
    'iii_cnn7_et_lg': (oof_cnn + oof_et + oof_lg) / 3.0,
}
sel = {k: roc_auc_score(ytr, v) for k, v in cand.items()}
best = max(sel, key=sel.get)
print('selection OOF aucs', {k: round(v, 4) for k, v in sel.items()}, '-> chose', best, flush=True)

# --- final: full-train fits, single test score ---
p_cnn = cnn_ens(train, test, SEED + 9000)
p_et = fit_et(train, ytr).predict_proba(dpc_arr(test))[:, 1]
p_lg = fit_lg(train, ytr).predict_proba(aac_arr(test))[:, 1]
final = {'i_cnn7': p_cnn, 'ii_cnn7_et': 0.5 * p_cnn + 0.5 * p_et,
         'iii_cnn7_et_lg': (p_cnn + p_et + p_lg) / 3.0}[best]
np.savez('/tmp/iter2_test.npz', p_cnn=p_cnn, p_et=p_et, p_lg=p_lg, final=final, yte=yte)

auc = roc_auc_score(yte, final); mcc = matthews_corrcoef(yte, final > 0.5)
rng = np.random.default_rng(SEED); n = len(yte); boot = []
for _ in range(10000):
    idx = rng.integers(0, n, n)
    if len(set(yte[idx])) < 2: continue
    boot.append(roc_auc_score(yte[idx], final[idx]))
ci = np.percentile(np.array(boot) - 0.83, [2.5, 97.5])
out = dict(model='iter2-' + best, dataset='anticp2_main OFFICIAL validation (n=%d)' % n,
           addendum='PREREG_ADDENDUM_2026-09-27.md iteration-2',
           selection_oof_aucs={k: float(v) for k, v in sel.items()},
           test_auc=float(auc), test_mcc=float(mcc),
           published_auc=0.83, published_mcc=0.51,
           bootstrap_diff_auc_minus_published_ci95=[float(ci[0]), float(ci[1])],
           beat=bool(auc > 0.83 and mcc > 0.51 and ci[0] > 0),
           elapsed_s=round(time.time() - t0, 1))
json.dump(out, open('results/stack_iter2.json', 'w'), indent=1)
print(json.dumps(out, indent=1), flush=True)
print('DONE', flush=True)
