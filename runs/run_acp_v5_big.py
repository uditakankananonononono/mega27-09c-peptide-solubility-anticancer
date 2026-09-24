"""ACP v5: wider/deeper CNNv2 (ch=128, 4 blocks, dil 1,2,4,8), 70 epochs,
cosine LR, label smoothing 0.05, plus a-priori equal-weight ensemble with
ET-DPC (weight fixed BEFORE seeing test)."""
import sys, json, time
sys.path.insert(0, 'src')
import numpy as np
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.metrics import roc_auc_score, accuracy_score, matthews_corrcoef, confusion_matrix
from pepx.datasets import augmented_acp_train, load_anticp2
from pepx.encoders import dpc
from pepx.models import PeptideCNNv2
from pepx.trainer import train_single_task

t0 = time.time()
train = augmented_acp_train('main')
test = [r for r in load_anticp2('main') if r.meta == 'test']
yte = np.array([r.label for r in test])

Xtr = np.stack([dpc(r.sequence) for r in train])
Xte = np.stack([dpc(r.sequence) for r in test])
ytr = np.array([r.label for r in train])
et = ExtraTreesClassifier(n_estimators=400, random_state=7, n_jobs=2).fit(Xtr, ytr)
p_et = et.predict_proba(Xte)[:, 1]
print('ET-DPC AUC', round(roc_auc_score(yte, p_et), 4), flush=True)

import torch.nn as nn
model = PeptideCNNv2(emb_dim=64, ch=128, blocks=4, dropout=0.25)
res, p_cnn = train_single_task(model, train, test, 'anticp2_main', 'CNNv2-big', epochs=70, wide=True, lr=1e-3, patience=15)
print(res, flush=True)

p_ens = 0.5 * p_et + 0.5 * p_cnn  # fixed a priori
pred = (p_ens >= 0.5).astype(int)
tn, fp, fn, tp = confusion_matrix(yte, pred, labels=[0, 1]).ravel()
out = dict(model='v5: ET-DPC + CNNv2-big equal-weight (a priori)', dataset='anticp2_main',
           cnnv2_big_auc=float(roc_auc_score(yte, p_cnn)), et_auc=float(roc_auc_score(yte, p_et)),
           accuracy=float(accuracy_score(yte, pred)), sensitivity=float(tp/(tp+fn)),
           specificity=float(tn/(tn+fp)), mcc=float(matthews_corrcoef(yte, pred)),
           auc=float(roc_auc_score(yte, p_ens)),
           published_anticp2=dict(acc=0.7399, mcc=0.48, auc=0.83))
print(out, flush=True)
json.dump(out, open('results/acp_v5.json', 'w'), indent=2)
print('elapsed', round(time.time()-t0, 1), 's', flush=True)
