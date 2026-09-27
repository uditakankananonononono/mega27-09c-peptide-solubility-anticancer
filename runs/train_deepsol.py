"""DeepSol solubility arm per PREREG_ADDENDUM_2026-09-27_SOLARM.md (locked).
Laptop build: in-RAM float32 features (saved to work/*.npy for resume), per-batch tensorization. Same locked
protocol: ET-DPC + CNNv2 (max_len 512, val early stop) + a priori mean;
single test evaluation per model; 10k bootstrap vs acc 0.77 / MCC 0.55.
"""
import sys, json, time
sys.path.insert(0, '.')
t0 = time.time()
import numpy as np
import torch
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.metrics import accuracy_score, matthews_corrcoef, roc_auc_score
from pepx.models import PeptideCNNv2
from pepx.encoders import dpc, aac, descriptor_array
from pepx.trainer import set_seed
from pepx.alphabet import AA_STANDARD

import os
os.makedirs('work', exist_ok=True)

SEED = 2709
set_seed(SEED)
PAD = len(AA_STANDARD); _AA = {a: i for i, a in enumerate(AA_STANDARD)}
MAXLEN = 512

def load(split):
    seqs = [l.strip() for l in open(f'data/{split}_src')]
    ys = np.array([int(l.strip()) for l in open(f'data/{split}_tgt')])
    return seqs, ys

tr_s, ytr = load('train'); va_s, yva = load('val'); te_s, yte = load('test')
print('sizes', len(tr_s), len(va_s), len(te_s), flush=True)

def feats_dpc(seqs, path):
    mm = np.zeros((len(seqs), 400), dtype='float32')
    for i in range(0, len(seqs), 2000):
        mm[i:i+2000] = np.stack([dpc(s) for s in seqs[i:i+2000]]).astype('float32')
    np.save(path, mm); return mm

def idx_map(seqs, path):
    mm = np.full((len(seqs), MAXLEN), PAD, dtype='int16')
    for i, s in enumerate(seqs):
        ids = [_AA.get(a, PAD) for a in s[:MAXLEN]]
        if ids: mm[i, :len(ids)] = ids
    np.save(path, mm); return mm

def wide_map(seqs, path):
    mm = np.zeros((len(seqs), 438), dtype='float32')
    for i in range(0, len(seqs), 2000):
        rows = []
        for s in seqs[i:i+2000]:
            rows.append(np.concatenate([descriptor_array(s), aac(s), dpc(s)]))
        mm[i:i+2000] = np.stack(rows).astype('float32')
    np.save(path, mm); return mm

import os
if os.path.exists('work/dsol_p_et_te.npy'):
    print('resume: ET stage already done, skipping', flush=True)
else:
    print('featurizing dpc...', flush=True)
    Xtr = feats_dpc(tr_s, 'work/dsol_xtr.npy')
    Xva_d = feats_dpc(va_s, 'work/dsol_xva_d.npy')
    Xte_d = feats_dpc(te_s, 'work/dsol_xte_d.npy')

    print('ET fit...', flush=True)
    et = ExtraTreesClassifier(n_estimators=400, random_state=SEED, n_jobs=2)
    et.fit(Xtr, ytr)
    p_et_va = et.predict_proba(Xva_d)[:, 1]
    p_et_te = et.predict_proba(Xte_d)[:, 1]
    print('ET val auc', round(roc_auc_score(yva, p_et_va), 4), flush=True)
    np.save('work/dsol_p_et_te.npy', p_et_te)
    np.save('work/dsol_p_et_va.npy', p_et_va)
    del Xtr, Xva_d, Xte_d, et
p_et_va = np.load('work/dsol_p_et_va.npy') if os.path.exists('work/dsol_p_et_va.npy') else None
p_et_te = np.load('work/dsol_p_et_te.npy')

if all(os.path.exists(p) for p in ['work/dsol_itr.npy','work/dsol_wtr.npy','work/dsol_iva.npy','work/dsol_wva.npy','work/dsol_ite.npy','work/dsol_wte.npy']):
    print('resume: wide+idx feature files exist, reusing', flush=True)
    Itr = np.load('work/dsol_itr.npy')
    Wtr = np.load('work/dsol_wtr.npy')
    Iva = np.load('work/dsol_iva.npy')
    Wva = np.load('work/dsol_wva.npy')
    Ite = np.load('work/dsol_ite.npy')
    Wte = np.load('work/dsol_wte.npy')
else:
    print('wide+idx featurizing...', flush=True)
    Itr = idx_map(tr_s, 'work/dsol_itr.npy')
    Wtr = wide_map(tr_s, 'work/dsol_wtr.npy')
    Iva = idx_map(va_s, 'work/dsol_iva.npy')
    Wva = wide_map(va_s, 'work/dsol_wva.npy')
    Ite = idx_map(te_s, 'work/dsol_ite.npy')
    Wte = wide_map(te_s, 'work/dsol_wte.npy')
mu = Wtr[:].mean(0, keepdims=True); sd = Wtr[:].std(0, keepdims=True).clip(1e-6)
np.save('work/dsol_mu.npy', mu); np.save('work/dsol_sd.npy', sd)

Hp = dict(emb_dim=32, ch=96, blocks=2, dropout=0.2780616350232075)
LR = 0.0017286476459924243
model = PeptideCNNv2(**Hp)
opt = torch.optim.Adam(model.parameters(), lr=LR)
lf = torch.nn.BCEWithLogitsLoss()
ytr_t = torch.tensor(ytr, dtype=torch.float32)
n = len(ytr_t)

def batch(i, d):
    xi = torch.from_numpy(np.asarray(i, dtype=np.int64))
    xd = torch.from_numpy((np.asarray(d, dtype=np.float32) - mu) / sd)
    return xi, xd

best_val, best_state, bad = -1.0, None, 0
start_ep = 0
if os.path.exists('work/dsol_full_ckpt.pt'):
    ck = torch.load('work/dsol_full_ckpt.pt')
    model.load_state_dict(ck['model']); opt.load_state_dict(ck['opt'])
    best_val, best_state, bad = ck['best_val'], ck['best_state'], ck['bad']
    start_ep = ck['epoch'] + 1
    print('resume: full checkpoint, next epoch', start_ep, 'best', round(best_val, 4), flush=True)
for ep in range(start_ep, 30):
    model.train(); perm = np.random.permutation(n)
    for i in range(0, n, 256):
        b = perm[i:i+256]
        xi, xd = batch(Itr[b], Wtr[b])
        opt.zero_grad()
        lf(model(xi, xd), ytr_t[torch.from_numpy(b)]).backward(); opt.step()
    model.eval()
    pv = []
    with torch.no_grad():
        for i in range(0, len(yva), 512):
            xi, xd = batch(Iva[i:i+512], Wva[i:i+512])
            pv.append(torch.sigmoid(model(xi, xd)))
    pv = torch.cat(pv).numpy()
    va_auc = roc_auc_score(yva, pv)
    print('epoch', ep, 'val auc', round(va_auc, 4), round(time.time() - t0), 's', flush=True)
    if va_auc > best_val:
        best_val, bad = va_auc, 0
        best_state = {k: v.clone() for k, v in model.state_dict().items()}
    else:
        bad += 1
    torch.save(dict(model=model.state_dict(), opt=opt.state_dict(), epoch=ep,
                    best_val=best_val, best_state=best_state, bad=bad),
               'work/dsol_full_ckpt.pt')
    if bad >= 3: print('early stop', flush=True); break
model.load_state_dict(best_state); model.eval()
pt = []
with torch.no_grad():
    for i in range(0, len(yte), 512):
        xi, xd = batch(Ite[i:i+512], Wte[i:i+512])
        pt.append(torch.sigmoid(model(xi, xd)))
p_cnn_te = torch.cat(pt).numpy()
np.save('work/dsol_p_cnn_te.npy', p_cnn_te)

def metrics(y, p):
    return dict(acc=float(accuracy_score(y, p > 0.5)), mcc=float(matthews_corrcoef(y, p > 0.5)),
                auc=float(roc_auc_score(y, p)))

rng = np.random.default_rng(SEED); m = len(yte)
def boot(fn, p, pub):
    vals = []
    for _ in range(10000):
        idx = rng.integers(0, m, m)
        vals.append(fn(yte[idx], p[idx]))
    return [float(x) for x in np.percentile(np.array(vals) - pub, [2.5, 97.5])]

res = {}
for name, p in [('ET-DPC', p_et_te), ('CNNv2', p_cnn_te), ('MEAN(ET+CNN)', (p_et_te + p_cnn_te) / 2.0)]:
    mt = metrics(yte, p)
    acc_ci = boot(lambda a, b: accuracy_score(a, b > 0.5), p, 0.77)
    mcc_ci = boot(lambda a, b: matthews_corrcoef(a, b > 0.5), p, 0.55)
    res[name] = dict(**mt, boot_acc_ci95=acc_ci, boot_mcc_ci95=mcc_ci,
                     beat=bool(mt['acc'] > 0.77 and mt['mcc'] > 0.55 and acc_ci[0] > 0 and mcc_ci[0] > 0))
    print(name, res[name], flush=True)

out = dict(addendum='PREREG_ADDENDUM_2026-09-27_SOLARM.md',
           partition='DeepSol official (Zenodo 10.5281/zenodo.1162886), test n=%d' % m,
           published=dict(acc=0.77, mcc=0.55),
           val_auc_et=float(roc_auc_score(yva, p_et_va)), val_auc_cnn=float(best_val),
           results=res, elapsed_s=round(time.time() - t0, 1))
json.dump(out, open('deepsol_arm.json', 'w'), indent=1)
print(json.dumps(out, indent=1), flush=True)
print('DONE', flush=True)
