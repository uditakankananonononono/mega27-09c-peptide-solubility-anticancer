"""WIP - NOT YET RUN. Figure batch for paper sections 15-16.
Tools: umap-learn, shap, logomaker, upsetplot, seaborn.
Each produces a committed figure + a small JSON of run provenance."""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(style="whitegrid")
from pepx.datasets import load_anticp2, load_cancerppd, load_esol
from pepx.encoders import descriptor_array
prov = {}

# 1. UMAP of descriptor space: accession panel + AntiCP2 sample + candidates
import umap
rows = json.load(open("results/accession_panels.json"))
KEYS = ["mw","charge","gravy","boman","aromaticity","instability"]
X_acc = np.array([[r[k] for k in KEYS] for r in rows])
labels = ["proteome"]*len(X_acc)
acp = load_anticp2("main")
import random; random.seed(7)
sub = random.sample(acp, 300)
from pepx.descriptors import net_charge, gravy, boman_index, aromaticity, instability_index, molecular_weight
def vec(s): return [molecular_weight(s), net_charge(s), gravy(s), boman_index(s), aromaticity(s), instability_index(s)]
X_acp = np.array([vec(d.sequence) for d in sub]); labels += ["anticp2_pos" if d.label else "anticp2_neg" for d in sub]
cands = ["FEKEAKKIEIKRH","KLAKLAKLAKLAK"]
X_c = np.array([vec(c) for c in cands]); labels += ["candidate"]*2
X = np.vstack([X_acc, X_acp, X_c])
mu, sd = X.mean(0), X.std(0)+1e-9
emb = umap.UMAP(n_neighbors=25, min_dist=0.4, random_state=7).fit_transform((X-mu)/sd)
plt.figure(figsize=(7,5))
for lab, mk in [("proteome","."),("anticp2_pos","."),("anticp2_neg","."),("candidate","*")]:
    m = np.array(labels)==lab
    plt.scatter(emb[m,0], emb[m,1], s=90 if lab=="candidate" else 8, marker=mk, label=lab, alpha=0.7)
plt.legend(); plt.title("UMAP of 6-descriptor space (z-scored)")
plt.tight_layout(); plt.savefig("paper/fig5_umap.png", dpi=160); plt.close()
prov["umap"] = {"n": int(len(X)), "keys": KEYS}

# 2. SHAP on RF solubility model
import shap
from sklearn.ensemble import RandomForestClassifier
es = load_esol()
Xe = np.array([descriptor_array(d.sequence) for d in es]); ye = np.array([d.label for d in es])
rf = RandomForestClassifier(n_estimators=200, random_state=7, n_jobs=2).fit(Xe, ye)
ex = shap.TreeExplainer(rf)
sv = ex.shap_values(Xe[:400])
sv1 = sv[1] if isinstance(sv, list) else sv[:,:,1] if sv.ndim==3 else sv
names = ["length","mw","net_charge","pI","gravy","eisenberg","hopp_woods","cf_helix","cf_sheet",
         "aggregation_zz","aromaticity","aliphatic_index","instability","boman",
         "frac_cationic","frac_anionic","frac_hydrophobic","frac_polar"]
imp = np.abs(sv1).mean(0)
order = np.argsort(imp)[::-1][:12]
plt.figure(figsize=(7,4))
plt.barh([names[i] for i in order][::-1], imp[order][::-1])
plt.title("SHAP mean |value| - RF solubility (eSOL)")
plt.tight_layout(); plt.savefig("paper/fig6_shap.png", dpi=160); plt.close()
prov["shap"] = {"top": [names[i] for i in order[:5]]}

# 3. logomaker: ACP positives vs negatives (positions 1-15)
import logomaker
pos = [d.sequence[:15] for d in acp if d.label==1 and len(d.sequence)>=15][:400]
neg = [d.sequence[:15] for d in acp if d.label==0 and len(d.sequence)>=15][:400]
def counts(rows):
    M = np.zeros((15,20))
    AA="ACDEFGHIKLMNPQRSTVWY"
    for r in rows:
        for i,a in enumerate(r): M[i, AA.index(a)] += 1
    import pandas as pd
    return pd.DataFrame(M, columns=list(AA))
fig, axes = plt.subplots(2,1, figsize=(8,5), sharex=True)
logomaker.Logo(counts(pos), ax=axes[0]); axes[0].set_title("ACP positives (n=400)")
logomaker.Logo(counts(neg), ax=axes[1]); axes[1].set_title("non-ACP (n=400)")
plt.tight_layout(); plt.savefig("paper/fig7_logo.png", dpi=160); plt.close()
prov["logomaker"] = {"pos": len(pos), "neg": len(neg)}

# 4. UpSet: sequence overlap across APD3 / AntiCP2 pos / CancerPPD / DBAASP
from upsetplot import from_memberships, plot as upset_plot
from pepx.datasets import load_apd3
sets = {"APD3": set(s[:20] for s in load_apd3() if len(s)>=20),
        "AntiCP2_pos": set(d.sequence[:20] for d in acp if d.label==1 and len(d.sequence)>=20),
        "CancerPPD": set(s[:20] for s in load_cancerppd() if len(s)>=20),
        "DBAASP": set(r["sequence"][:20] for r in json.load(open("data/external/dbaasp_records.json")) if len(r["sequence"])>=20)}
from collections import Counter
cnt = Counter(tuple(k for k,v in sets.items() if s0 in v) for s0 in set().union(*sets.values()))
data = from_memberships(list(cnt.keys()), data=list(cnt.values()))
plt.figure(figsize=(7,4)); upset_plot(data, sort_by="cardinality", show_counts=True)
plt.savefig("paper/fig8_upset.png", dpi=160, bbox_inches="tight"); plt.close()
prov["upset"] = {k: len(v) for k,v in sets.items()}

json.dump(prov, open("results/tool_figures.json","w"), indent=1)
print("figures done:", list(prov))
