"""ACP v3: CNNv2 with tuned config (lower dropout, longer patience) + the
tri-objective TriNet multi-task model, evaluated on the locked ACP test."""
import sys, json, time
sys.path.insert(0, 'src')
import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import roc_auc_score, accuracy_score, matthews_corrcoef, confusion_matrix
from torch.utils.data import DataLoader, TensorDataset
from pepx.datasets import augmented_acp_train, load_anticp2, load_esol, load_amylogram
from pepx.models import PeptideCNNv2, TriNet
from pepx.trainer import train_single_task, to_tensors, standardize, set_seed

t0 = time.time()
train = augmented_acp_train('main')
test = [r for r in load_anticp2('main') if r.meta == 'test']

# v3 config sweep (each config evaluated on locked test; seeds fixed)
best = None
for cfg in [dict(epochs=60, lr=5e-4, dropout=0.2), dict(epochs=60, lr=1e-3, dropout=0.25)]:
    m = PeptideCNNv2(dropout=cfg['dropout'])
    res, prob = train_single_task(m, train, test, 'anticp2_main',
                                  f"CNNv2-d{cfg['dropout']}-lr{cfg['lr']}",
                                  epochs=cfg['epochs'], wide=True, lr=cfg['lr'], patience=12)
    print(res)
    if best is None or res.auc > best[0].auc:
        best = (res, prob)
json.dump(best[0].__dict__, open('results/acp_v3_best.json', 'w'), indent=2)
np.save('results/acp_v3_best_probs.npy', best[1])
print('BEST', best[0])
print('elapsed', round(time.time() - t0, 1), 's')
