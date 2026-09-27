"""C8 random-weight CNN falsifier (addendum 2026-09-27, verdict section 7:
"Random-weight CNN falsifier not run"). Same AntiCP2 main test set, same
encoding/standardization as the trained CNNv2, but weights stay at seeded
random init (seed 2709, no training). Expectation if training carries the
signal: near-chance metrics. Reported as-is either way."""
import sys, json
sys.path.insert(0, 'src')
import numpy as np, torch
from sklearn.metrics import accuracy_score, matthews_corrcoef, roc_auc_score
from pepx.models import PeptideCNNv2
from pepx.datasets import augmented_acp_train, load_anticp2
from pepx.trainer import to_tensors, standardize, set_seed

set_seed(2709)
train = augmented_acp_train('main')
test = [r for r in load_anticp2('main') if r.meta == 'test']
Xtr_i, Xtr_d, ytr = to_tensors(train, 60, wide=True)
Xte_i, Xte_d, yte = to_tensors(test, 60, wide=True)
_, Xte_d = standardize(Xtr_d, Xte_d)

model = PeptideCNNv2()   # random init, never trained
model.eval()
with torch.no_grad():
    p = torch.sigmoid(model(Xte_i, Xte_d)).numpy()
y = yte.numpy()
res = {
 'protocol': 'random-weight falsifier; identical encoding/standardization/test split as cnnv2_aug_anticp2_main; seed 2709; no training',
 'n_test': int(len(y)),
 'random_weight': dict(acc=float(accuracy_score(y, p > 0.5)),
                       mcc=float(matthews_corrcoef(y, p > 0.5)),
                       auc=float(roc_auc_score(y, p))),
 'trained_reference': json.load(open('results/cnnv2_aug_anticp2_main.json')),
}
json.dump(res, open('results/random_weight_falsifier.json', 'w'), indent=1)
print(json.dumps(res['random_weight'], indent=1))
