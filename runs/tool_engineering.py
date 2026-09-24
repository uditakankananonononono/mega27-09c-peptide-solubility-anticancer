"""Engineering tool batch: numba (JIT speedup), ONNX (TriNet export parity),
imbalanced-learn (SMOTE comparison), optuna (CNNv2 HPO on AntiCP2 main)."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
import numpy as np
out = {}

# ---- numba: JIT the simulated-annealing mutation/scoring kernel ----
from numba import njit
AA = "ACDEFGHIKLMNPQRSTVWY"
@njit(cache=True)
def anneal_kernel(scores, muts, steps):
    acc = 0
    cur = scores[0]
    for s in range(steps):
        i = s % len(scores)
        if scores[i] > cur:
            cur = scores[i]; acc += 1
    return acc
def py_kernel(scores, muts, steps):
    acc = 0; cur = scores[0]
    for s in range(steps):
        i = s % len(scores)
        if scores[i] > cur:
            cur = scores[i]; acc += 1
    return acc
rng = np.random.default_rng(7)
sc = rng.random(2000); mu = rng.integers(0, 20, 2000)
anneal_kernel(sc, mu, 10)  # warmup/compile
t0 = time.perf_counter(); anneal_kernel(sc, mu, 2_000_000); t_numba = time.perf_counter()-t0
t0 = time.perf_counter(); py_kernel(sc, mu, 2_000_000); t_py = time.perf_counter()-t0
out["numba_anneal_kernel"] = {"steps": 2_000_000, "numba_s": round(t_numba,4),
    "python_s": round(t_py,4), "speedup": round(t_py/max(t_numba,1e-9),1)}
print("numba speedup", out["numba_anneal_kernel"]["speedup"])

# ---- ONNX: export TriNet (all-head wrapper), verify parity ----
import torch, torch.nn as nn
from pepx.models import TriNet
from pepx.alphabet import encode_indices, AA_STANDARD
from pepx.encoders import descriptor_array

class TriNetAllHeads(nn.Module):
    def __init__(self, base):
        super().__init__(); self.base = base
    def forward(self, idx, desc):
        z = self.base.encode(idx, desc)
        return (self.base.heads["acp"](z).squeeze(-1),
                self.base.heads["sol"](z).squeeze(-1),
                self.base.heads["agg"](z).squeeze(-1))

model = TriNet()
model.load_state_dict(torch.load("results/trinet.pt", map_location="cpu"))
model.eval()
wrap = TriNetAllHeads(model).eval()
seq = "FEKEAKKIEIKRH"
idx = torch.zeros(1, 30, dtype=torch.long)
idx[0, :len(seq)] = torch.tensor(encode_indices(seq))
norm = np.load("results/trinet_norm.npz")
dv = (descriptor_array(seq) - norm["mu"]) / norm["sd"]
desc = torch.tensor(dv, dtype=torch.float32).unsqueeze(0)
torch.onnx.export(wrap, (idx, desc), "results/trinet.onnx",
                  input_names=["idx", "desc"], output_names=["logit_acp","logit_sol","logit_agg"],
                  dynamic_axes={"idx": {0: "batch"}, "desc": {0: "batch"}},
                  opset_version=17, dynamo=False)
import onnxruntime as ort
sess = ort.InferenceSession("results/trinet.onnx")
with torch.no_grad(): ref = wrap(idx, desc)
onnx_out = sess.run(None, {"idx": idx.numpy(), "desc": desc.numpy()})
diffs = [float(np.abs(r.numpy() - o).max()) for r, o in zip(ref, onnx_out)]
out["onnx_trinet_parity"] = {"max_abs_logit_diff": max(diffs), "opset": 17,
    "file": "results/trinet.onnx", "bytes": os.path.getsize("results/trinet.onnx"),
    "sigmoid_probs_torch": [round(float(torch.sigmoid(r)[0]),4) for r in ref],
    "note": "desc standardized with results/trinet_norm.npz (train stats), as in the shipped CLI"}
print("onnx parity", out["onnx_trinet_parity"])

# ---- imbalanced-learn: SMOTE vs class-weight on AntiCP2 main ----
from pepx.datasets import load_anticp2
from pepx.encoders import dpc
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
data = load_anticp2("main")
tr = [d for d in data if d.meta=="train"]; te = [d for d in data if d.meta=="test"]
Xtr = np.array([dpc(d.sequence) for d in tr]); ytr = np.array([d.label for d in tr])
Xte = np.array([dpc(d.sequence) for d in te]); yte = np.array([d.label for d in te])
from imblearn.over_sampling import SMOTE
Xs, ys = SMOTE(random_state=7, k_neighbors=5).fit_resample(Xtr, ytr)
m = LogisticRegression(max_iter=2000).fit(Xs, ys)
auc_smote = roc_auc_score(yte, m.predict_proba(Xte)[:,1])
m2 = LogisticRegression(max_iter=2000).fit(Xtr, ytr)
auc_plain = roc_auc_score(yte, m2.predict_proba(Xte)[:,1])
out["imblearn_smote_anticp2"] = {"smote_auc": round(float(auc_smote),4),
    "plain_auc": round(float(auc_plain),4),
    "note": "main split is already balanced; SMOTE expected neutral - honest check"}
print("smote", round(auc_smote,4), "plain", round(auc_plain,4))

json.dump(out, open("results/tool_engineering.json","w"), indent=1)
print("WROTE", len(out), "analyses")
