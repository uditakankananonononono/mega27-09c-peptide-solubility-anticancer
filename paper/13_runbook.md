# Reproduction runbook

Every number in this paper regenerates with the commands below, in dependency order, on the reference container (2 CPU / 1.9 GB / no GPU; Python 3.10; `pip install -r requirements.txt` equivalents noted per section). Total wall time: under 3 hours. Random seeds are fixed at 7 everywhere; where a result is single-seed, Section 18 (Limitations) says so and sizes the variance band.

## 0. Test gate

```
python3 -m pytest          # 29/29 green, ~7 s, hermetic (no network)
```

Covers alphabet/validation contracts, FASTA round-trips, descriptor identities (pI monotonicity, Ikai identity, water-subtraction identity), the modlAMP mass lock, the full 20-residue pyteomics mass lock, encoder shapes, CLI end-to-end, and 900 hypothesis property cases.

## 1. Data acquisition

```
python3 runs/fetch_uniprot_accessions.py   # 128 UniProt entries -> data/uniprot_accessions/
python3 runs/fetch_external_dbs.py         # ChEMBL (68 records), PDB (8 entries), attempt logs for CAMP/Hemolytik
python3 runs/fetch_dbaasp_novelty.py       # 169 DBAASP records + novelty screen
```

AntiCP 2.0, CancerPPD, eSOL, APD3, AmyloGram/pep424 raw files: downloaded once from the URLs in `data/MANIFEST.md` (each entry carries the live source URL and the local path). No raw file in this repository was hand-edited; validation is code, not curation.

## 2. Descriptor validation

```
python3 -c "from pepx.descriptors import descriptor_vector; print(descriptor_vector('KLAKLAKLAKLAK'))"
python3 runs/accession_panels.py           # 118-entry proteome panel -> results/accession_panels.json
```

Cross-validation hierarchy: modlAMP (first validator), peptides.py (second), pyteomics (third - the one that caught the P/Q/V/W/Y scramble, Section 4), propy3 (CTD family, rank correlation). All outputs in `results/tool_analyses.json` and `results/descriptor_crossvalidation.json`.

## 3. Benchmarks

```
python3 runs/run_baselines.py              # classical grid, both AntiCP2 splits
python3 runs/run_cnnv2_acp.py              # CNNv2+augmentation on main
python3 runs/run_ensemble_acp.py           # 4-arm ensemble (negative)
python3 runs/run_solubility.py             # eSOL, locked seed-7 split
python3 runs/run_aggregation.py            # AmyloGram GNN
python3 runs/run_pep424_v3.py              # THE head-to-head: transfer + CV5
python3 runs/run_optuna_acp.py             # 12-trial HPO on main (~2 min)
```

v1/v2 pep424 scripts are kept (invalid joins, documented) - running them reproduces the *wrong* numbers on purpose; the paper cites them as the negative trail.

## 4. Tool analyses and figures

```
python3 runs/tool_analyses.py              # xgboost/lightgbm/statsmodels/networkx/propy3/pyteomics/biotite
python3 runs/tool_engineering.py           # numba/ONNX/SMOTE
python3 runs/tool_figures.py               # fig5-fig8 (UMAP/SHAP/logos/UpSet)
python3 runs/pdb_chembl_analysis.py        # structure metrics + potency table
python3 runs/api_smoke.py                  # FastAPI TestClient smoke
```

## 5. Paper

```
cd paper && pandoc 00_front.md 01_math.md 02_data.md 03_methods_results.md \
  04_discovery_tool.md 05_related.md 06_appendix.md 07_proteomics.md \
  08_figures_analysis.md 09_engineering.md 10_limitations.md \
  11_math_extended.md 12_datasets.md 13_runbook.md \
  -o paper.pdf --pdf-engine=pdflatex && pdfinfo paper.pdf
```

Page counts in all external reports come from `pdfinfo` on the built artifact - never from source-length estimates.
