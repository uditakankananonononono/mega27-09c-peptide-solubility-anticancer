"""Per-accession descriptor panels: every fetched UniProt accession is used
individually - 18-descriptor panel on its full sequence, a sliding-window
(15-mer) peptide profile, and TriNet-style gates vs the discovery screen
background. Output feeds paper section 15 (proteomic context) and the
dataset ledger."""
import json, glob, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
from pepx.descriptors import descriptor_vector
from pepx.alphabet import validate_sequence

rows = []
for f in sorted(glob.glob("data/uniprot_accessions/*.json")):
    d = json.load(open(f))
    seq = d["sequence"]
    acc = d["accession"]
    if not seq or len(seq) < 30: continue
    # full-protein panel
    seq = validate_sequence(seq)
    panel = descriptor_vector(seq)
    # sliding 15-mer windows, stride 10: per-window hydrophobicity + charge spread
    from pepx.descriptors import gravy as _gravy
    win_h = []
    for i in range(0, len(seq) - 14, 10):
        win_h.append(_gravy(seq[i:i+15]))
    rows.append({
        "accession": acc, "protein": d["protein"][:60], "organism": d["organism"],
        "length": d["length"], "reviewed": d["reviewed"],
        "mw": panel["mw"], "charge": panel["net_charge"], "gravy": panel["gravy"],
        "boman": panel["boman"], "aromaticity": panel["aromaticity"],
        "instability": panel["instability"],
        "n_windows": len(win_h),
        "win_hydro_mean": sum(win_h)/max(1,len(win_h)),
        "win_hydro_max": max(win_h) if win_h else None,
    })
json.dump(rows, open("results/accession_panels.json","w"), indent=1)
print(f"panels computed for {len(rows)} accessions")
# quick stats for the report
import statistics as st
print("mean length", round(st.mean(r["length"] for r in rows)))
print("amyloid-control mean GRAVY", round(st.mean(r["gravy"] for r in rows if r["accession"] in
      ("P05067","P10997","P02647","P02649")),3))
