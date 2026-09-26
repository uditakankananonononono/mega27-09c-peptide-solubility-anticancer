"""Iteration-1 stacked ensemble per PREREG_ADDENDUM_2026-09-27 (locked before scoring).
Base: (a) CNNv2 HPO config, (b) ET-DPC(400), (c) AAC+train-mined-motif logistic.
Meta: logistic on 3 OOF probability columns. Single official-test evaluation +
10,000-replicate bootstrap vs published AUROC 0.83 / MCC 0.51.
"""
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
set_seed(SEED)
recs = load_anticp2('main')
train = [r for r in recs if r.meta == 'train']
test = [r for r in recs if r.meta == 'test']
ytr = np.array([r.label for r in train])
yte = np.array([r.label for r in test])
print('train', len(train), 'test', len(test), flush=True)

MAXLEN = 50
Hp = dict(emb_dim=32, ch=96, blocks=2, dropout=0.2780616350232075)
LR = 0.0017286476459924243
EPOCHS = 18

def cnn_fit_predict(tr_recs, ev_recs, seed):
    set_seed(seed)
    model = PeptideCNNv2(**Hp)
    Xi, Xd, y = to_tensors(tr_recs, MAXLEN, wide=True)
    Ei, Ed, _ = to_tensors(ev_recs, MAXLEN, wide=True)
    Xd, Ed = standardize(Xd, Ed)
    opt = torch.optim.Adam(model.parameters(), lr=LR)
    lossf = torch.nn.BCEWithLogitsLoss()
    n = len(y)
    for ep in range(EPOCHS):
        model.train()
        perm = torch.randperm(n)
        for i in range(0, n, 64):
            b = perm[i:i+64]
            opt.zero_grad()
            loss = lossf(model(Xi[b], Xd[b]), y[b])
            loss.backward(); opt.step()
    model.eval()
    with torch.no_grad():
        p = torch.sigmoid(model(Ei, Ed)).numpy()
    return p

def dpc_arr(rs): return np.stack([dpc(r.sequence) for r in rs])
def aac_arr(rs): return np.stack([aac(r.sequence) for r in rs])

def mine_motifs(rs, ys, k_range=(4, 5, 6), top=20):
    from collections import Counter
    pos_c, neg_c = Counter(), Counter()
    for r, y in zip(rs, ys):
        s = r.sequence
        for k in k_range:
            for i in range(len(s) - k + 1):
                (pos_c if y == 1 else neg_c)[s[i:i+k]] += 1
    scored = []
    for m, c in pos_c.items():
        if c >= 5:
            prec = c / (c + neg_c.get(m, 0) + 1)
            scored.append((prec, c, m))
    scored.sort(reverse=True)
    return [m for _, _, m in scored[:top]]

def motif_feats(rs, motifs):
    X = np.zeros((len(rs), len(motifs)))
    for i, r in enumerate(rs):
        for j, m in enumerate(motifs):
            if m in r.sequence: X[i, j] = 1
    return X

oof = np.zeros((len(train), 3))
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
for fold, (ti, vi) in enumerate(skf.split(train, ytr)):
    tr_r = [train[i] for i in ti]; va_r = [train[i] for i in vi]
    y_ti = ytr[ti]
    print('fold', fold, 'cnn...', flush=True)
    oof[vi, 0] = cnn_fit_predict(tr_r, va_r, SEED + fold)
    et = ExtraTreesClassifier(n_estimators=400, random_state=SEED, n_jobs=2)
    et.fit(dpc_arr(tr_r), y_ti)
    oof[vi, 1] = et.predict_proba(dpc_arr(va_r))[:, 1]
    motifs = mine_motifs(tr_r, y_ti)
    Xm = np.hstack([aac_arr(tr_r), motif_feats(tr_r, motifs)])
    Xv = np.hstack([aac_arr(va_r), motif_feats(va_r, motifs)])
    lg = LogisticRegression(max_iter=2000)
    lg.fit(Xm, y_ti)
    oof[vi, 2] = lg.predict_proba(Xv)[:, 1]
    print('fold', fold, 'aucs', [round(roc_auc_score(ytr[vi], oof[vi, j]), 4) for j in range(3)], flush=True)

np.savez('/tmp/stack_iter1_oof.npz', oof=oof, ytr=ytr)
meta = LogisticRegression(max_iter=2000)
meta.fit(oof, ytr)
print('OOF stack auc', roc_auc_score(ytr, meta.predict_proba(oof)[:, 1]), flush=True)

# full-train refits -> single test evaluation
p_test = np.zeros((len(test), 3))
p_test[:, 0] = cnn_fit_predict(train, test, SEED + 100)
et = ExtraTreesClassifier(n_estimators=400, random_state=SEED, n_jobs=2)
et.fit(dpc_arr(train), ytr)
p_test[:, 1] = et.predict_proba(dpc_arr(test))[:, 1]
motifs = mine_motifs(train, ytr)
lg = LogisticRegression(max_iter=2000)
lg.fit(np.hstack([aac_arr(train), motif_feats(train, motifs)]), ytr)
p_test[:, 2] = lg.predict_proba(np.hstack([aac_arr(test), motif_feats(test, motifs)]))[:, 1]
s_test = meta.predict_proba(p_test)[:, 1]
np.savez('/tmp/stack_iter1_test.npz', p_test=p_test, s_test=s_test, yte=yte)

auc = roc_auc_score(yte, s_test)
mcc = matthews_corrcoef(yte, s_test > 0.5)
rng = np.random.default_rng(SEED)
n = len(yte)
boot_auc = []
for _ in range(10000):
    idx = rng.integers(0, n, n)
    if len(set(yte[idx])) < 2: continue
    boot_auc.append(roc_auc_score(yte[idx], s_test[idx]))
boot_auc = np.array(boot_auc)
ci = np.percentile(boot_auc - 0.83, [2.5, 97.5])
out = dict(model='iter1-stack(CNNv2-hpo + ET-DPC + motif-logreg, OOF meta-logreg)',
           dataset='anticp2_main OFFICIAL validation (n=%d)' % n,
           addendum='PREREG_ADDENDUM_2026-09-27.md',
           oof_stack_auc=float(roc_auc_score(ytr, meta.predict_proba(oof)[:, 1])),
           test_auc=float(auc), test_mcc=float(mcc),
           published_auc=0.83, published_mcc=0.51,
           bootstrap_diff_auc_minus_published_ci95=[float(ci[0]), float(ci[1])],
           beat=bool(auc > 0.83 and mcc > 0.51 and ci[0] > 0),
           elapsed_s=round(time.time() - t0, 1))
json.dump(out, open('results/stack_iter1.json', 'w'), indent=1)
print(json.dumps(out, indent=1), flush=True)
print('DONE', flush=True)
