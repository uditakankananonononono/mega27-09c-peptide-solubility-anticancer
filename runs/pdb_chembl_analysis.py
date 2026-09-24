"""PDB (biotite) + ChEMBL analysis: structure context for known ACPs and
external bioactivity evidence table. Feeds paper section 17."""
import json, urllib.request, os
import numpy as np

def get(url, timeout=60):
    req = urllib.request.Request(url, headers={"User-Agent":"pepx-research/0.2 (academic)"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()

out = {}
# --- biotite structure metrics on ACP structures (first NMR model) ---
import biotite.structure.io.pdb as bpdb
import biotite.structure as bs
for pid in ["2K6O","2MAG"]:
    try:
        path = f"data/external/{pid}.pdb"
        if not os.path.exists(path):
            open(path,"wb").write(get(f"https://files.rcsb.org/download/{pid}.pdb"))
        arr = bpdb.get_structure(bpdb.PDBFile.read(path), model=1)
        ca = arr[arr.atom_name == "CA"]
        coord = ca.coord
        rg = float(np.sqrt(((coord - coord.mean(0))**2).sum(1).mean()))
        d = np.linalg.norm(coord[1:] - coord[:-1], axis=1)
        out[pid] = {"n_res": int(len(ca)), "radius_of_gyration_A": round(rg,2),
                    "mean_ca_ca_A": round(float(d.mean()),2),
                    "chain_ids": sorted(set(ca.chain_id.tolist()))}
    except Exception as e:
        out[pid] = {"error": str(e)[:100]}
print(out)

# --- ChEMBL bioactivity table ---
d = json.load(open("data/external/chembl_acp_bioactivity.json"))
rows = {}
for a in d["activities"]:
    mol = a["mol"]
    if a.get("value") is None: continue
    try: v = float(a["value"])
    except (TypeError, ValueError): continue
    rows.setdefault(mol, []).append({"type": a["type"], "value": v, "units": a["units"], "target": a["target"]})
table = {m: {"n_records": len(v),
             "potency_range": [min(x["value"] for x in v), max(x["value"] for x in v)],
             "units_seen": sorted(set(x["units"] or "" for x in v)),
             "targets": sorted(set(x["target"] or "" for x in v))[:6]}
         for m, v in rows.items()}
out["chembl_table"] = table
json.dump(out, open("results/pdb_chembl_analysis.json","w"), indent=1)
print("molecules with activity:", len(table))
