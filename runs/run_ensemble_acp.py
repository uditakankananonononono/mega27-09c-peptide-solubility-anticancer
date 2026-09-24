"""Train all model arms on AntiCP2 main, tune ensemble weights on the
validation carve only, evaluate the ensemble ONCE on the locked test split."""
import sys, json, time
sys.path.insert(0, 'src')
import numpy as np
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.metrics import roc_auc_score, accuracy_score, matthews_corrcoef, confusion_matrix
from pepx.datasets import augmented_acp_train, load_anticp2
from pepx.encoders import dpc
from pepx.models import PeptideCNN, PeptideCNNv2, PeptideGNN
from pepx.trainer import train_single_task, to_tensors, standardize

t0 = time.time()
train = augmented_acp_train('main')
test = [r for r in load_anticp2('main') if r.meta == 'test']

# --- ET-DPC arm (classical champion) ---
Xtr = np.stack([dpc(r.sequence) for r in train])
Xte = np.stack([dpc(r.sequence) for r in test])
ytr = np.array([r.label for r in train])
yte = np.array([r.label for r in test])
et = ExtraTreesClassifier(n_estimators=400, random_state=7, n_jobs=2).fit(Xtr, ytr)
p_et = et.predict_proba(Xte)[:, 1]

# --- deep arms ---
res_cnn, p_cnn = train_single_task(PeptideCNNv2(), train, test, 'anticp2_main', 'CNNv2', epochs=40, wide=True)
res_gnn, p_gnn = train_single_task(PeptideGNN(), train, test, 'anticp2_main', 'GNN', epochs=30)
res_c1, p_c1 = train_single_task(PeptideCNN(), train, test, 'anticp2_main', 'CNN', epochs=30)

arms = {'ET-DPC': p_et, 'CNNv2': p_cnn, 'GNN': p_gnn, 'CNN': p_c1}
for name, p in arms.items():
    pred = (p >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(yte, pred, labels=[0, 1]).ravel()
    print(f"{name:8s} acc={accuracy_score(yte, pred):.4f} mcc={matthews_corrcoef(yte, pred):.4f} auc={roc_auc_score(yte, p):.4f}")

# --- probability ensemble (equal + greedy AUC weights chosen on TEST would be
# cheating; use equal weights and report) ---
p_ens = np.mean(np.stack(list(arms.values())), axis=0)
pred = (p_ens >= 0.5).astype(int)
tn, fp, fn, tp = confusion_matrix(yte, pred, labels=[0, 1]).ravel()
out = dict(model='ENSEMBLE(ET+CNNv2+GNN+CNN)', dataset='anticp2_main',
           accuracy=float(accuracy_score(yte, pred)),
           sensitivity=float(tp / (tp + fn)), specificity=float(tn / (tn + fp)),
           mcc=float(matthews_corrcoef(yte, pred)), auc=float(roc_auc_score(yte, p_ens)))
print(out)
json.dump({'ensemble': out,
           'arms': {k: {'auc': float(roc_auc_score(yte, v))} for k, v in arms.items()}},
          open('results/ensemble_anticp2_main.json', 'w'), indent=2)
np.savez('results/probs_anticp2_main.npz', yte=yte, **arms)
print('elapsed', round(time.time() - t0, 1), 's')
