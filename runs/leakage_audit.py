"""Train/test leakage audit for AntiCP2 main: k-mer containment similarity
between every test peptide and its nearest train neighbor."""
import sys, json
sys.path.insert(0, 'src')
from pepx.datasets import load_anticp2
data = load_anticp2('main')
tr = [d.sequence for d in data if d.meta == 'train']
te = [d.sequence for d in data if d.meta == 'test']
def kmers(s, k=4): return {s[i:i+k] for i in range(len(s)-k+1)} if len(s) >= k else {s}
tr_k = [(s, kmers(s)) for s in tr]
sims = []
exact = 0
for t in te:
    tk = kmers(t)
    best = 0.0
    for s, sk in tr_k:
        u = len(tk | sk)
        if u == 0: continue
        j = len(tk & sk) / u
        if j > best: best = j
    if t in tr: exact += 1
    sims.append(best)
import numpy as np
sims = np.array(sims)
out = {"n_test": len(te), "n_train": len(tr), "k": 4,
       "exact_duplicates": exact,
       "max_jaccard_mean": round(float(sims.mean()),4),
       "frac_test_above_0.5": round(float((sims > 0.5).mean()),4),
       "frac_test_above_0.7": round(float((sims > 0.7).mean()),4),
       "frac_test_above_0.9": round(float((sims > 0.9).mean()),4)}
json.dump(out, open("results/leakage_audit_anticp2.json","w"), indent=1)
print(out)
