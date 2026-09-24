"""eSOL threshold sensitivity: RF-DESC AUROC at solubility cutoffs 20-50%."""
import sys, json
sys.path.insert(0, 'src')
import numpy as np
from pepx.datasets import load_esol
from pepx.encoders import descriptor_array
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
out = []
for thr in (20, 30, 40, 50):
    es = load_esol(threshold_pct=thr)
    X = np.array([descriptor_array(d.sequence) for d in es]); y = np.array([d.label for d in es])
    i_tr, i_te = train_test_split(np.arange(len(y)), test_size=0.15, stratify=y, random_state=7)
    m = RandomForestClassifier(n_estimators=200, random_state=7, n_jobs=2).fit(X[i_tr], y[i_tr])
    auc = roc_auc_score(y[i_te], m.predict_proba(X[i_te])[:,1])
    out.append({"threshold_pct": thr, "n": int(len(y)), "pos_frac": round(float(y.mean()),3),
                "rf_desc_auc": round(float(auc),4)})
    print(out[-1])
json.dump(out, open("results/esol_threshold_sensitivity.json","w"), indent=1)
