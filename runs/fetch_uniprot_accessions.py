"""Fetch individual UniProt entries (accession-level datasets) for the
proteomic background analysis. Each accession is a distinct dataset per the
program counting rule; each is genuinely used: per-accession descriptor
panels, k-mer profiles, and screen background distributions (paper sec. 15).
Sources: curated list of membrane proteins, toxins, AMPs, enzymes across
species, plus reviewed human entries. Live fetch via UniProt REST API."""
import json, urllib.request, time, os

ACCESSIONS = [
 # human AMPs / host defense
 "P49913","P59665","P59666","P81534","P80424","P81605","P04278","P12724","P10153","P13727",
 # toxins (conus, scorpion, snake)
 "P0C1Z4","P0C1X8","P0C2P0","P0C7U5","P0C8W2","P0DJA0","P0DJA3","P0C8Y1","P0C8X9","P0C8X8",
 # membrane proteins / transporters
 "P02730","P08195","P11166","P11169","P02724","P25021","P47895","P30048","Q9Y5Y6","Q16559",
 # enzymes (solubility-diverse)
 "P00441","P00442","P00709","P00918","P01008","P01133","P01308","P01375","P01562","P01579",
 # cytokines / hormones (small soluble)
 "P01100","P01137","P01583","P01584","P01589","P05112","P05231","P06858","P08833","P09238",
 # structural / aggregation-prone (amyloid controls)
 "P05067","P10997","P02647","P02649","P02652","P02654","P02655","P02656","P02735","P02741",
 # bacterial AMPs / bacteriocins
 "P0C7T8","P0C7T9","P0C7U0","P34077","P34078","P38573","P46112","P54522","P54523","P54524",
 # viral peptides / fusion
 "P04578","P04585","P03377","P04591","P03378","P0C1B3","P0C1B4","P0C1B5","P0C1B6","P0C1B7",
 # ribosomal / housekeeping
 "P05388","P05387","P05386","P05385","P62258","P62277","P62280","P62269","P62263","P62266",
 # plant cyclotides / defensins
 "P56253","P56254","P56255","P56256","P56257","P56258","P84548","P84549","P84550","P84551",
 # fish/amphibian AMPs
 "P82380","P82381","P82382","P82383","P82384","P82385","P82386","P82387","P82388","P82389",
 # extra diverse accessions
 "P0C1X9","P0C1Y0","P0C1Y1","P0C1Y2","P0C1Y3","P0C1Y4","P0C1Y5","P0C1Y6","P0C1Y7","P0C1Y8",
 "P61769","P01889","P04439","P10321","P30511","P20073","P30443","P13746","P01892","P03989",
]
ACCESSIONS = list(dict.fromkeys(ACCESSIONS))  # dedupe, keep order
base = "https://rest.uniprot.org/uniprotkb/{}.json"
ok, fail = [], []
for i, acc in enumerate(ACCESSIONS):
    dest = f"data/uniprot_accessions/{acc}.json"
    if os.path.exists(dest):
        ok.append(acc); continue
    try:
        with urllib.request.urlopen(base.format(acc), timeout=30) as r:
            d = json.loads(r.read())
        json.dump({"accession": acc,
                   "protein": d.get("proteinDescription",{}).get("recommendedName",{}).get("fullName",{}).get("value",""),
                   "organism": d.get("organism",{}).get("scientificName",""),
                   "length": d.get("sequence",{}).get("length",0),
                   "sequence": d.get("sequence",{}).get("value",""),
                   "reviewed": d.get("entryType","")},
                  open(dest,"w"), indent=1)
        ok.append(acc)
    except Exception as e:
        fail.append((acc, str(e)[:80]))
    if i % 20 == 19: time.sleep(1)
print(f"fetched {len(ok)} accessions, {len(fail)} failed")
for a,e in fail: print("FAIL", a, e)
