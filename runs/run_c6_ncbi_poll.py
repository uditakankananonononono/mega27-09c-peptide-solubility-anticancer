#!/usr/bin/env python3
"""C6: poll NCBI BLAST RIDs and parse top-hit identity/coverage.

Prereg 5d3103d novelty rule: NOVEL = no hit with >=90% identity over >=80%
of query length. Resumable: results written incrementally to
results/c6_partial/ncbi_results.json. NCBI etiquette: one request per RID,
>=30s between requests, per prereg.
"""
import json, os, re, time, urllib.request, urllib.parse

RIDS = "results/c6_partial/ncbi_rids.json"
OUT = "results/c6_partial/ncbi_results.json"
WAIT_S = 30.0
MAX_ROUNDS = 40  # safety bound; each round re-polls unfinished RIDs

def fetch(rid):
    url = ("https://blast.ncbi.nlm.nih.gov/Blast.cgi?" +
           urllib.parse.urlencode({"CMD": "Get", "RID": rid,
                                   "FORMAT_TYPE": "Text"}))
    req = urllib.request.Request(url, headers={"User-Agent": "c6-novelty-audit/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode("utf-8", "replace")

def parse_top_hit(text, qlen):
    """Return dict of best alignment: pident, coverage, accession, desc.
    Parse the first alignment block: 'Identities = a/b (p%)' and the Query
    coordinate span across its lines."""
    m = re.search(r"Sequences producing significant alignments", text)
    if not m:
        return None  # no hits
    # first alignment block
    blocks = re.split(r"\n>+", text)
    if len(blocks) < 2:
        return None
    blk = blocks[1]
    desc = blk.splitlines()[0].strip()[:200]
    acc = desc.split()[0] if desc else ""
    idents = re.findall(r"Identities\s*=\s*(\d+)/(\d+)\s*\((\d+)%\)", blk)
    if not idents:
        return None
    # take the first (highest-scoring) HSP
    a, b, p = idents[0]
    pident = int(p) / 100.0
    # query span across this HSP's Query lines (up to next HSP 'Score =' or end)
    hsp = blk.split("Score =", 2)
    seg = hsp[1] if len(hsp) > 1 else blk
    qcoords = re.findall(r"Query\s+(\d+)\s+\S+\s+(\d+)", seg)
    if qcoords:
        qstart = min(int(s) for s, _ in qcoords)
        qend = max(int(e) for _, e in qcoords)
        cov = (qend - qstart + 1) / qlen
    else:
        cov = None
    return {"pident": pident, "coverage": cov, "accession": acc, "desc": desc}

def main():
    rids = json.load(open(RIDS))
    results = json.load(open(OUT)) if os.path.exists(OUT) else {}
    pending = [k for k in rids if k not in results or results[k].get("status") == "waiting"]
    print(f"{len(rids)} RIDs, {len(pending)} pending", flush=True)
    for rnd in range(MAX_ROUNDS):
        if not pending:
            break
        still = []
        for key in pending:
            rid = rids[key]["rid"]
            try:
                txt = fetch(rid)
            except Exception as e:
                print(f"{key} fetch error {e}", flush=True)
                still.append(key)
                time.sleep(WAIT_S)
                continue
            if "Status=WAITING" in txt:
                results[key] = {"status": "waiting", "rid": rid}
                still.append(key)
            elif "Status=UNKNOWN" in txt or "Status=FAILED" in txt:
                results[key] = {"status": "error", "rid": rid,
                                "note": "RID unknown/failed at NCBI; needs resubmit"}
            else:
                qlen = len(key.split(":", 1)[1])
                hit = parse_top_hit(txt, qlen)
                rec = {"status": "done", "rid": rid}
                if hit is None:
                    rec["top_hit"] = None
                    rec["novel"] = True
                else:
                    hit["not_novel"] = (hit["pident"] >= 0.90 and
                                        (hit["coverage"] or 0) >= 0.80)
                    rec["top_hit"] = hit
                    rec["novel"] = not hit["not_novel"]
                results[key] = rec
                print(f"{key} done novel={rec['novel']}", flush=True)
            json.dump(results, open(OUT, "w"), indent=1)
            time.sleep(WAIT_S)
        pending = still
        if pending:
            print(f"round {rnd}: {len(pending)} still waiting", flush=True)
    json.dump(results, open(OUT, "w"), indent=1)
    print("POLL DONE", flush=True)

if __name__ == "__main__":
    main()
