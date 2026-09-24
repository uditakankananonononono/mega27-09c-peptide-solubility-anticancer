"""Fetch individual DBAASP peptide records by ID (accession-level datasets)
and run the candidate novelty screen against them (falsifiability criterion i)."""
import json, urllib.request, time, ssl

def get(url, timeout=40):
    req = urllib.request.Request(url, headers={"User-Agent":"pepx-research/0.2 (academic)"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()

# individual record fetches: monomer ids 36..236 (DBAASPR ids resolved per record)
records, fails = [], 0
for pid in range(36, 237):
    try:
        d = json.loads(get(f"https://dbaasp.org/peptides/{pid}?format=json"))
        rec = d.get("data", d)
        if isinstance(rec, list): rec = rec[0] if rec else {}
        seq = rec.get("sequence") or ""
        if seq:
            records.append({"id": pid, "dbaaspId": rec.get("dbaaspId"),
                            "name": rec.get("name"), "sequence": seq,
                            "length": rec.get("sequenceLength")})
        else:
            fails += 1
    except Exception:
        fails += 1
    if pid % 40 == 0: time.sleep(0.5)
json.dump(records, open("data/external/dbaasp_records.json","w"), indent=1)
print(f"dbaasp individual records: {len(records)} fetched, {fails} skipped/failed")

# novelty screen: exact + substring match of candidates against all external pools
cands = ["FEKEAKKIEIKRH"]
try:
    d = json.load(open("results/discovery_candidates_scaled.json"))
    cands += [c["sequence"] for c in d[:20]]
except Exception: pass
pool = [r["sequence"] for r in records]
try:
    d0 = json.load(open("data/external/dbaasp_sample.json"))
    for rec in d0.get("data", []):
        for m in rec.get("monomers", []) or []:
            if m.get("sequence"): pool.append(m["sequence"])
except Exception: pass
hits = []
for c in cands:
    for p in pool:
        if c == p or (len(c) >= 8 and c in p) or (len(p) >= 8 and p in c):
            hits.append({"candidate": c, "match": p[:40]})
json.dump({"candidates_screened": len(cands), "pool_size": len(pool), "hits": hits},
          open("results/dbaasp_novelty_screen.json","w"), indent=1)
print(f"novelty: {len(cands)} candidates vs pool {len(pool)} -> {len(hits)} hits")
