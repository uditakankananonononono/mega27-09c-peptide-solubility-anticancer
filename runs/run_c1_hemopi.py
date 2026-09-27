"""C1 hemolysis fifth task (PREREG_C1_HEMOLYSIS_2026-09-27.md, locked before
scoring). PeptideCNNv2 + train_single_task AntiCP2 recipe, no augmentation;
train on published main split, score once on validation. Falsifier:
label-shuffle model on HemoPI-1 (seed 13)."""
import sys, json, time
sys.path.insert(0, 'src')
import numpy as np
from pathlib import Path
from pepx.datasets import LabeledPeptide
from pepx.models import PeptideCNNv2
from pepx.trainer import train_single_task

RAW = Path('data/raw/hemopi')

def load_hemopi(ds):
    recs = []
    for part in ('main', 'validation'):
        for cls, lab in (('pos', 1), ('neg', 0)):
            f = RAW / f'hemopi{ds}_{part}_{cls}.fa'
            seq = None
            for line in open(f):
                line = line.strip()
                if line.startswith('>'):
                    if seq:
                        recs.append(LabeledPeptide(seq, lab, f'hemopi{ds}', part))
                    seq = ''
                else:
                    seq = (seq or '') + line
            if seq:
                recs.append(LabeledPeptide(seq, lab, f'hemopi{ds}', part))
    return recs

def bootstrap_mcc(labels, probs, seed=23, n=10000):
    from sklearn.metrics import matthews_corrcoef
    labels = np.asarray(labels); preds = (np.asarray(probs) >= 0.5).astype(int)
    rng = np.random.RandomState(seed)
    stats = []
    N = len(labels)
    for _ in range(n):
        idx = rng.randint(0, N, N)
        stats.append(matthews_corrcoef(labels[idx], preds[idx]))
    return float(np.percentile(stats, 2.5)), float(np.percentile(stats, 97.5))

def decision(mcc, lo, hi, bar):
    if mcc > bar and lo > bar:
        return 'BEAT'
    if lo <= bar <= hi:
        return 'MATCH'
    return 'MISS'

t0 = time.time()
out = {'prereg': 'docs/PREREG_C1_HEMOLYSIS_2026-09-27.md',
       'comparators': {'hemopi1': {'mcc': 0.93, 'sn': 96.4, 'sp': 99.1, 'acc': 96.4},
                       'hemopi2': {'mcc': 0.51, 'sn': 78.2, 'sp': 78.3, 'acc': 75.7}},
       'datasets': {}}
for ds, bar in ((1, 0.93), (2, 0.51), (3, None)):
    recs = load_hemopi(ds)
    train = [r for r in recs if r.meta == 'main']
    val = [r for r in recs if r.meta == 'validation']
    res, prob = train_single_task(PeptideCNNv2(), train, val, f'hemopi{ds}',
                                  'PeptideCNNv2', epochs=40, wide=True)
    labels = [r.label for r in val]
    lo, hi = bootstrap_mcc(labels, prob)
    d = {'mcc': res.mcc, 'auc': res.auc, 'accuracy': res.accuracy,
         'sensitivity': res.sensitivity, 'specificity': res.specificity,
         'mcc_ci95': [lo, hi], 'n_val': len(val),
         'epochs_run': res.epochs_run, 'best_val_auc': res.best_val_auc}
    if bar is not None:
        d['comparator_mcc'] = bar
        d['decision'] = decision(res.mcc, lo, hi, bar)
    out['datasets'][f'hemopi{ds}'] = d
    np.save(f'results/c1_hemopi{ds}_probs.npy', prob)
    print(f'hemopi{ds}: mcc {res.mcc:.4f} ci [{lo:.4f},{hi:.4f}] auc {res.auc:.4f} {d.get("decision","")}', flush=True)

# falsifier: label-shuffle on hemopi1 (seed 13), same recipe
recs = load_hemopi(1)
train = [r for r in recs if r.meta == 'main']
val = [r for r in recs if r.meta == 'validation']
rng = np.random.RandomState(13)
labels = np.array([r.label for r in train])
perm = rng.permutation(labels)
train_shuf = [LabeledPeptide(r.sequence, int(perm[i]), r.source, r.meta)
              for i, r in enumerate(train)]
res_f, prob_f = train_single_task(PeptideCNNv2(), train_shuf, val, 'hemopi1',
                                  'PeptideCNNv2-shuffled', epochs=40, wide=True)
out['falsifier'] = {'shuffle_seed': 13, 'val_mcc': res_f.mcc,
                    'abs_mcc_lt_0.15': bool(abs(res_f.mcc) < 0.15)}
print(f"falsifier: shuffled-label val MCC {res_f.mcc:.4f}", flush=True)
json.dump(out, open('results/c1_hemopi.json', 'w'), indent=2)
print('elapsed', round(time.time() - t0, 1), 's')
