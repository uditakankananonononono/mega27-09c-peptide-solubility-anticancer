"""C4 CamSol-style per-residue profile (PREREG_C4_CAMSOL_PROFILE_2026-09-27.md).
4-weight ridge fit on eSOL seed-7 train only; Spearman + AUC@30% on test;
shuffled-label falsifier (seed 13)."""
import sys, json
sys.path.insert(0, 'src')
import numpy as np
from pepx.datasets import load_esol
from pepx.descriptors import KD_HYDROPATHY, CF_HELIX, CF_SHEET

AA = "ACDEFGHIKLMNPQRSTVWY"
CHARGE = dict(zip(AA, [0, -1, -1, -1, 0, 0, 0.1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0]))

def zscale(d):
    v = np.array([d[a] for a in AA], dtype=float)
    v = (v - v.mean()) / v.std()
    return dict(zip(AA, v))

ZH, ZA, ZB, ZC = zscale(KD_HYDROPATHY), zscale(CF_HELIX), zscale(CF_SHEET), zscale(CHARGE)

def profile(seq, w=(1, 1, 1, 1)):
    """Per-residue linear combination, then 7-residue window mean (edge-padded)."""
    s = np.array([w[0]*ZH[a] + w[1]*ZC[a] + w[2]*ZA[a] + w[3]*ZB[a] for a in seq])
    k = 7; pad = k // 2
    sp = np.concatenate([np.full(pad, s[0]), s, np.full(pad, s[-1])])
    return np.convolve(sp, np.ones(k)/k, mode='valid')

def features(seq):
    """Window-mean of each scale (linearity makes the fit exact for the form)."""
    out = []
    for Z in (ZH, ZC, ZA, ZB):
        s = np.array([Z[a] for a in seq]); k = 7; pad = k // 2
        sp = np.concatenate([np.full(pad, s[0]), s, np.full(pad, s[-1])])
        out.append(np.convolve(sp, np.ones(k)/k, mode='valid').mean())
    return out

def sol_pct(r):
    return float(r.meta.split('=')[1])

recs = load_esol()
idx = json.load(open('results/esol_split_seed7.json'))
tr = [recs[i] for i in idx['train_idx']]
te = [recs[i] for i in idx['test_idx']]
Xtr = np.array([features(r.sequence) for r in tr])
ytr = np.array([sol_pct(r) for r in tr])
Xte = np.array([features(r.sequence) for r in te])
yte = np.array([sol_pct(r) for r in te])

def fit(X, y, ridge=1.0):
    Xb = np.hstack([X, np.ones((len(X), 1))])
    A = Xb.T @ Xb; A[np.diag_indices_from(A)] += ridge
    return np.linalg.solve(A, Xb.T @ y)

def spearman(a, b):
    from scipy.stats import spearmanr
    return float(spearmanr(a, b).statistic)

beta = fit(Xtr, ytr)
w, b0 = beta[:4], beta[4]
score_te = Xte @ w + b0
rho = spearman(score_te, yte)
from sklearn.metrics import roc_auc_score
lab_te = (yte >= 30.0).astype(int)
auc = roc_auc_score(lab_te, score_te)
rng = np.random.RandomState(23)
boots = []
for _ in range(10000):
    ii = rng.randint(0, len(yte), len(yte))
    if len(set(lab_te[ii])) < 2:
        continue
    boots.append(roc_auc_score(lab_te[ii], score_te[ii]))
lo, hi = float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))
bar = 0.7915
decision = 'BEAT' if (auc > bar and lo > bar) else ('MATCH' if lo <= bar <= hi else 'MISS')
# falsifier
beta_f = fit(Xtr, np.random.RandomState(13).permutation(ytr))
rho_f = spearman(Xte @ beta_f[:4] + beta_f[4], yte)
examples = {}
for name, seq in (('c7_top_candidate', 'KLFKKICRQAEKDKDS'),
                  ('esol_test_first', te[0].sequence)):
    p = profile(seq, w)
    examples[name] = {'length': len(seq), 'profile_mean': float(p.mean()),
                      'profile_min': float(p.min()), 'profile_max': float(p.max()),
                      'profile': [round(float(x), 4) for x in p[:60]]}
out = {'prereg': 'docs/PREREG_C4_CAMSOL_PROFILE_2026-09-27.md',
       'disclosure': 'CamSol-style re-implementation; Table S3 weights inaccessible (paywalled SI); weights fit in-house on eSOL seed-7 train split ONLY; pattern/gatekeeper corrections not implemented (v1 scope cut, locked).',
       'fitted_weights': {'a_H_KD': float(w[0]), 'a_C_charge': float(w[1]),
                          'a_alpha_CF': float(w[2]), 'a_beta_CF': float(w[3]), 'intercept': float(b0)},
       'test': {'n': len(te), 'spearman': rho, 'auc_30pct': float(auc),
                'auc_ci95': [lo, hi], 'comparator_rf_desc_auc': bar, 'decision': decision},
       'falsifier': {'seed': 13, 'shuffled_spearman': rho_f, 'abs_lt_0.1': bool(abs(rho_f) < 0.1)},
       'examples': examples}
json.dump(out, open('results/c4_camsol_profile.json', 'w'), indent=2)
print(json.dumps({k: out[k] for k in ('fitted_weights', 'test', 'falsifier')}, indent=1))
