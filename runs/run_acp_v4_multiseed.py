"""ACP v4: multi-seed CNNv2 probability-averaged ensemble on AntiCP2 main.
Seeds chosen a priori (7, 17, 27); single evaluation on the locked test."""
import sys, json, time
sys.path.insert(0, 'src')
import numpy as np
from sklearn.metrics import roc_auc_score, accuracy_score, matthews_corrcoef, confusion_matrix
from pepx.datasets import augmented_acp_train, load_anticp2
from pepx.models import PeptideCNNv2
from pepx.trainer import train_single_task

t0 = time.time()
train = augmented_acp_train('main')
test = [r for r in load_anticp2('main') if r.meta == 'test']
yte = np.array([r.label for r in test])

probs = []
for seed in (7, 17, 27):
    res, p = train_single_task(PeptideCNNv2(dropout=0.25), train, test,
                               'anticp2_main', f'CNNv2-s{seed}', epochs=50,
                               wide=True, lr=7e-4, patience=12, seed=seed)
    print(res, flush=True)
    probs.append(p)
p_ens = np.mean(np.stack(probs), axis=0)
pred = (p_ens >= 0.5).astype(int)
tn, fp, fn, tp = confusion_matrix(yte, pred, labels=[0, 1]).ravel()
out = dict(model='CNNv2-seedens(7,17,27)', dataset='anticp2_main',
           accuracy=float(accuracy_score(yte, pred)),
           sensitivity=float(tp / (tp + fn)), specificity=float(tn / (tn + fp)),
           mcc=float(matthews_corrcoef(yte, pred)), auc=float(roc_auc_score(yte, p_ens)))
print(out, flush=True)
json.dump(out, open('results/acp_v4_seedens.json', 'w'), indent=2)
np.save('results/acp_v4_seedens_probs.npy', p_ens)
print('elapsed', round(time.time() - t0, 1), 's', flush=True)
