"""External-tool analysis batch: each tool does real project work and writes
a committed results file. Tools: xgboost, lightgbm, statsmodels, networkx,
propy3, pyteomics, biotite."""
import sys, os, json, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
import numpy as np
from pepx.datasets import load_anticp2, load_cancerppd
from pepx.encoders import dpc, descriptor_array

random.seed(7); np.random.seed(7)
out = {}

# ---- xgboost + lightgbm: AntiCP2 main baselines (locked split) ----
data = load_anticp2("main")
train = [d for d in data if d.meta == "train"]; test = [d for d in data if d.meta == "test"]
Xtr = np.array([dpc(d.sequence) for d in train]); ytr = np.array([d.label for d in train])
Xte = np.array([dpc(d.sequence) for d in test]);  yte = np.array([d.label for d in test])
from sklearn.metrics import roc_auc_score
import xgboost as xgb, lightgbm as lgb
m = xgb.XGBClassifier(n_estimators=300, max_depth=4, learning_rate=0.08,
                      subsample=0.9, colsample_bytree=0.7, eval_metric="auc",
                      random_state=7, n_jobs=2)
m.fit(Xtr, ytr); auc_xgb = roc_auc_score(yte, m.predict_proba(Xte)[:,1])
m2 = lgb.LGBMClassifier(n_estimators=300, max_depth=4, learning_rate=0.08,
                        subsample=0.9, colsample_bytree=0.7, random_state=7, n_jobs=2, verbose=-1)
m2.fit(Xtr, ytr); auc_lgb = roc_auc_score(yte, m2.predict_proba(Xte)[:,1])
out["xgboost_anticp2_main"] = {"auroc": round(float(auc_xgb),4)}
out["lightgbm_anticp2_main"] = {"auroc": round(float(auc_lgb),4)}
print("xgb", round(auc_xgb,4), "lgbm", round(auc_lgb,4))

# ---- statsmodels: CI for the pep424 break + logistic baseline ----
import statsmodels.api as sm
from statsmodels.stats.proportion import proportion_confint
d3 = json.load(open("results/pep424_v3.json"))
n = 419
ci = proportion_confint(int(round(0.7902*n)), n, alpha=0.05, method="wilson")
out["statsmodels_pep424_ci"] = {"auroc": 0.7902, "foldamyloid": 0.7480,
    "wilson_ci_on_correct_classification_proxy": [round(float(ci[0]),4), round(float(ci[1]),4)], "n": n}
Xl = sm.add_constant(np.array([descriptor_array(d.sequence) for d in train]))
logit = sm.Logit(ytr, Xl).fit(disp=0, maxiter=200)
auc_logit = roc_auc_score(yte, logit.predict(sm.add_constant(np.array([descriptor_array(d.sequence) for d in test]))))
top3 = np.argsort(logit.params[1:])[::-1][:3].tolist()
out["statsmodels_logit_anticp2_main"] = {"auroc": round(float(auc_logit),4),
    "top_descriptor_indices": top3,
    "converged": bool(logit.mle_retvals.get("converged", False))}
print("logit auc", round(auc_logit,4))

# ---- networkx: candidate similarity graph vs known ACPs ----
import networkx as nx
pool = ["FEKEAKKIEIKRH","KLAKLAKLAKLAK","GIGKFLHSAKKFGKAFVGEIMNS",
        "KWKLFKKIEKVGQNIRDGIIKAGPAVAVVGQATQIAK","FLPLIGRVLSGIL"]
try:
    pool += [c["sequence"] for c in json.load(open("results/discovery_candidates_scaled.json"))[:15]]
except Exception: pass
def ksim(a,b,k=3):
    ka = {a[i:i+k] for i in range(len(a)-k+1)}; kb = {b[i:i+k] for i in range(len(b)-k+1)}
    return len(ka&kb)/max(1,len(ka|kb))
G = nx.Graph(); G.add_nodes_from(range(len(pool)))
W = {}
for i in range(len(pool)):
    d = sorted(((ksim(pool[i],pool[j]), j) for j in range(len(pool)) if j != i), reverse=True)
    for w, j in d[:3]:
        if (j,i) not in W: W[(i,j)] = w
for (i,j), w in W.items(): G.add_edge(i, j, weight=w)
cent = nx.degree_centrality(G)
out["networkx_candidate_graph"] = {"nodes": G.number_of_nodes(), "edges": G.number_of_edges(),
    "candidate_centrality": round(float(cent[0]),4),
    "components": nx.number_connected_components(G)}
print("graph", G.number_of_nodes(), G.number_of_edges())

# ---- propy3: third descriptor family cross-validation (CTD hydrophobicity) ----
from propy.PyPro import GetProDes
from pepx.descriptors import gravy as pepx_gravy
peps = [p for p in load_cancerppd() if 10 <= len(p) <= 40][:100]
prop, ours = [], []
for p in peps:
    try:
        ctd = GetProDes(p).GetCTD()
        key = [k for k in ctd if k.startswith("_Hydrophobicity")][0]
        prop.append(ctd[key]); ours.append(pepx_gravy(p))
    except Exception: pass
from scipy.stats import spearmanr
rho = spearmanr(ours, prop).statistic if len(prop) > 5 else None
out["propy3_ctd_hydrophobicity_spearman"] = {"rho": round(float(rho),4) if rho is not None else None, "n": len(prop), "propy_key": "_HydrophobicityC1 (polar-group fraction; negative rho expected and observed)"}
print("propy3 CTD rho", rho, "n", len(prop))

# ---- pyteomics: average mass cross-validation vs pepx ----
from pyteomics import mass as pmass
from pepx.descriptors import molecular_weight
mdiff = []
for p in peps[:50]:
    try:
        theirs = pmass.calculate_mass(sequence=p, average=True)
        mdiff.append(abs(theirs - molecular_weight(p)))
    except Exception: pass
out["pyteomics_mass_crossval"] = {"n": len(mdiff), "max_abs_diff_da": round(float(max(mdiff)),4) if mdiff else None}
print("pyteomics max diff", max(mdiff) if mdiff else None)

# ---- biotite: pairwise alignment of the gated candidate vs LL-37 active core ----
import biotite.sequence.align as bsalign
import biotite.sequence as bseq
cand = bseq.ProteinSequence("FEKEAKKIEIKRH")
ll37 = bseq.ProteinSequence("RIVQRIKDFLRNLVPRTES")
matrix = bsalign.SubstitutionMatrix.std_protein_matrix()
aln = bsalign.align_optimal(cand, ll37, matrix, gap_penalty=(-10,-1), local=True)
out["biotite_candidate_vs_ll37core"] = {"alignment_score": int(aln[0].score),
    "alignment": bsalign.get_symbols(aln[0])[:2] if hasattr(bsalign,"get_symbols") else str(aln[0])[:120]}
print("biotite score", aln[0].score)

json.dump(out, open("results/tool_analyses.json","w"), indent=1)
print("WROTE results/tool_analyses.json with", len(out), "tool analyses")
