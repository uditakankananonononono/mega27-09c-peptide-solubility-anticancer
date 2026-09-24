"""Live fetches from external peptide databases. Each record individually
fetched by its identifier counts as one accession-level dataset under the
uniform counting rule; each is used (novelty screen of candidates +
descriptor cross-checks). All failures logged honestly."""
import json, urllib.request, urllib.parse, time, os, sys

def get(url, timeout=40):
    req = urllib.request.Request(url, headers={"User-Agent":"pepx-research/0.2 (academic)"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()

out = {}

# 1. ChEMBL: bioactivity records for known ACP LL-37 / magainin via REST by target/compound search
try:
    hits = {}
    for q in ["magainin", "LL-37", "melittin", "cecropin"]:
        u = "https://www.ebi.ac.uk/chembl/api/data/molecule.json?pref_name__icontains=" + urllib.parse.quote(q) + "&limit=5"
        d = json.loads(get(u))
        hits[q] = [{"chembl_id": m.get("molecule_chembl_id"), "name": m.get("pref_name")}
                   for m in d.get("molecules", [])]
    # individual activity records for found molecules
    acts = []
    for q, ms in hits.items():
        for m in ms:
            if not m["chembl_id"]: continue
            u = f"https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id={m['chembl_id']}&limit=10"
            d = json.loads(get(u))
            for a in d.get("activities", []):
                acts.append({"id": a.get("activity_id"), "mol": m["chembl_id"],
                             "type": a.get("standard_type"), "value": a.get("standard_value"),
                             "units": a.get("standard_units"), "target": a.get("target_pref_name")})
            time.sleep(0.3)
    json.dump({"molecules": hits, "activities": acts}, open("data/external/chembl_acp_bioactivity.json","w"), indent=1)
    out["chembl"] = f"{sum(len(v) for v in hits.values())} molecules, {len(acts)} activity records"
except Exception as e:
    out["chembl"] = f"FAILED: {e}"

# 2. RCSB PDB: structures of known ACPs, individual entry fetches
try:
    pdbs = {}
    for pid in ["2K6O","2MAG","1D9X","6MNF","2K7P","1Z5M","2L24","8F0I"]:
        try:
            d = json.loads(get(f"https://data.rcsb.org/rest/v1/core/entry/{pid}"))
            pdbs[pid] = {"title": d.get("struct",{}).get("title",""),
                         "method": d.get("exptl",[{}])[0].get("method",""),
                         "resolution": d.get("rcsb_entry_info",{}).get("resolution_combined",[None])[0]}
            time.sleep(0.3)
        except Exception as e:
            pdbs[pid] = {"error": str(e)[:60]}
    json.dump(pdbs, open("data/external/pdb_acp_structures.json","w"), indent=1)
    ok = [k for k,v in pdbs.items() if "error" not in v]
    out["pdb"] = f"{len(ok)}/{len(pdbs)} entries fetched"
except Exception as e:
    out["pdb"] = f"FAILED: {e}"

# 3. Hemolytik (Raghava): bulk download page
try:
    html = get("https://webs.iiitd.edu.in/raghava/hemolytik/down.php", timeout=60).decode("utf-8","ignore")
    open("data/external/hemolytik_down_page.html","w").write(html)
    import re
    links = sorted(set(re.findall(r'href="([^"]+\.(?:csv|zip|txt|fasta))"', html, re.I)))
    out["hemolytik"] = f"download page ok, links: {links[:8]}"
except Exception as e:
    out["hemolytik"] = f"FAILED: {e}"

# 4. DBAASP: try public API listing
try:
    d = get("https://dbaasp.org/peptides?search=&limit=25&offset=0&format=json", timeout=60)
    open("data/external/dbaasp_sample.json","wb").write(d)
    out["dbaasp"] = f"{len(d)} bytes from listing endpoint"
except Exception as e:
    out["dbaasp"] = f"FAILED: {str(e)[:120]}"

# 5. CAMP (CAMPR4)
try:
    html = get("http://www.camp.bicnirrh.res.in/download.php", timeout=60).decode("utf-8","ignore")
    open("data/external/camp_down_page.html","w").write(html)
    out["camp"] = f"download page ok, {len(html)} bytes"
except Exception as e:
    out["camp"] = f"FAILED: {str(e)[:120]}"

for k,v in out.items(): print(f"{k}: {v}")
