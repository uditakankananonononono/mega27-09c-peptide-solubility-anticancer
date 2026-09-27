"""DeepSol solubility arm per PREREG_ADDENDUM_2026-09-27_SOLARM.md (locked).
Memory-lean: float32 memmap features, per-batch tensorization. Same locked
protocol: ET-DPC + CNNv2 (max_len 512, val early stop) + a priori mean;
single test evaluation per model; 10k bootstrap vs acc 0.77 / MCC 0.55.
"""
import sys, json, time
sys.path.insert(0, 'src')
t0 = time.time()
import numpy as np
import torch
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.metrics import accuracy_score, matthews_corrcoef, roc_auc_score
from pepx.models import PeptideCNNv2
from pepx.encoders import dpc, aac, descriptor_array
from pepx.trainer import set_seed
from pepx.alphabet import AA_STANDARD

SEED = 2709
set_seed(SEED)
PAD = len(AA_STANDARD); _AA = {a: i for i, a in enumerate(AA_STANDARD)}
MAXLEN = 512

def load(split):
    seqs = [l.strip() for l in open(f'data/raw/deepsol/{split}_src')]
    ys = np.array([int(l.strip()) for l in open(f'data/raw/deepsol/{split}_tgt')])
    return seqs, ys

tr_s, ytr = load('train'); va_s, yva = load('val'); te_s, yte = load('test')
print('sizes', len(tr_s), len(va_s), len(te_s), flush=True)

def feats_dpc(seqs, path):
    mm = np.memmap(path, dtype='float32', mode='w+', shape=(len(seqs), 400))
    for i in range(0, len(seqs), 2000):
        mm[i:i+2000] = np.stack([dpc(s) for s in seqs[i:i+2000]]).astype('float32')
    mm.flush(); return mm

def idx_map(seqs, path):
    mm = np.memmap(path, dtype='int16', mode='w+', shape=(len(seqs), MAXLEN))
    mm[:] = PAD
    for i, s in enumerate(seqs):
        ids = [_AA.get(a, PAD) for a in s[:MAXLEN]]
        if ids: mm[i, :len(ids)] = ids
    mm.flush(); return mm

def wide_map(seqs, path):
    mm = np.memmap(path, dtype='float32', mode='w+', shape=(len(seqs), 438))
    for i in range(0, len(seqs), 2000):
        rows = []
        for s in seqs[i:i+2000]:
            rows.append(np.concatenate([descriptor_array(s), aac(s), dpc(s)]))
        mm[i:i+2000] = np.stack(rows).astype('float32')
    mm.flush(); return mm

import os
if os.path.exists('/tmp/dsol_p_et_te.npy'):
    print('resume: ET stage already done, skipping', flush=True)
else:
    print('featurizing dpc...', flush=True)
    Xtr = feats_dpc(tr_s, '/tmp/dsol_xtr.mmap')
    Xva_d = feats_dpc(va_s, '/tmp/dsol_xva_d.mmap')
    Xte_d = feats_dpc(te_s, '/tmp/dsol_xte_d.mmap')

    print('ET fit...', flush=True)
    et = ExtraTreesClassifier(n_estimators=400, random_state=SEED, n_jobs=2)
    et.fit(Xtr, ytr)
    p_et_va = et.predict_proba(Xva_d)[:, 1]
    p_et_te = et.predict_proba(Xte_d)[:, 1]
    print('ET val auc', round(roc_auc_score(yva, p_et_va), 4), flush=True)
    np.save('/tmp/dsol_p_et_te.npy', p_et_te)
    np.save('/tmp/dsol_p_et_va.npy', p_et_va)
    del Xtr, Xva_d, Xte_d, et
p_et_va = np.load('/tmp/dsol_p_et_va.npy') if os.path.exists('/tmp/dsol_p_et_va.npy') else None
p_et_te = np.load('/tmp/dsol_p_et_te.npy')

if all(os.path.exists(p) for p in ['/tmp/dsol_itr.mmap','/tmp/dsol_wtr.mmap','/tmp/dsol_iva.mmap','/tmp/dsol_wva.mmap','/tmp/dsol_ite.mmap','/tmp/dsol_wte.mmap']):
    print('resume: wide+idx memmaps exist, reusing', flush=True)
    Itr = np.memmap('/tmp/dsol_itr.mmap', dtype='int16', mode='r', shape=(len(tr_s), MAXLEN))
    Wtr = np.memmap('/tmp/dsol_wtr.mmap', dtype='float32', mode='r', shape=(len(tr_s), 438))
    Iva = np.memmap('/tmp/dsol_iva.mmap', dtype='int16', mode='r', shape=(len(va_s), MAXLEN))
    Wva = np.memmap('/tmp/dsol_wva.mmap', dtype='float32', mode='r', shape=(len(va_s), 438))
    Ite = np.memmap('/tmp/dsol_ite.mmap', dtype='int16', mode='r', shape=(len(te_s), MAXLEN))
    Wte = np.memmap('/tmp/dsol_wte.mmap', dtype='float32', mode='r', shape=(len(te_s), 438))
else:
    print('wide+idx featurizing...', flush=True)
    Itr = idx_map(tr_s, '/tmp/dsol_itr.mmap')
    Wtr = wide_map(tr_s, '/tmp/dsol_wtr.mmap')
    Iva = idx_map(va_s, '/tmp/dsol_iva.mmap')
    Wva = wide_map(va_s, '/tmp/dsol_wva.mmap')
    Ite = idx_map(te_s, '/tmp/dsol_ite.mmap')
    Wte = wide_map(te_s, '/tmp/dsol_wte.mmap')
mu = Wtr[:].mean(0, keepdims=True); sd = Wtr[:].std(0, keepdims=True).clip(1e-6)
np.save('/tmp/dsol_mu.npy', mu); np.save('/tmp/dsol_sd.npy', sd)

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
if os.path.exists('/tmp/dsol_full_ckpt.pt'):
    ck = torch.load('/tmp/dsol_full_ckpt.pt')
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
               '/tmp/dsol_full_ckpt.pt')
    if bad >= 3: print('early stop', flush=True); break
model.load_state_dict(best_state); model.eval()
pt = []
with torch.no_grad():
    for i in range(0, len(yte), 512):
        xi, xd = batch(Ite[i:i+512], Wte[i:i+512])
        pt.append(torch.sigmoid(model(xi, xd)))
p_cnn_te = torch.cat(pt).numpy()
np.save('/tmp/dsol_p_cnn_te.npy', p_cnn_te)

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
json.dump(out, open('results/deepsol_arm.json', 'w'), indent=1)
print(json.dumps(out, indent=1), flush=True)
print('DONE', flush=True)
