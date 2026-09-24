# Appendix D. Per-configuration record

Every training run executed during this project, in execution order, with its committed result. "Locked test" metrics are single evaluations; "val" is the stratified carve used for early stopping.

| # | script | dataset | model | key metric | verdict |
|---|---|---|---|---|---|
| 1 | scripts_run_baselines.py | anticp2_main | SVM-AAC | 0.579 AUC | floor |
| 2 | scripts_run_baselines.py | anticp2_main | RF-DPC | 0.779 AUC | floor |
| 3 | scripts_run_baselines.py | anticp2_main | ET-DPC | 0.799 AUC | classical best |
| 4 | scripts_run_baselines.py | anticp2_alternate | ET-DPC | 0.931 AUC | below bar |
| 5 | run_cnnv2_acp.py | anticp2_main | CNNv2+APD3-aug | 0.781 AUC | augmentation neutral |
| 6 | run_ensemble_acp.py | anticp2_main | CNNv2 | 0.801 AUC | deep best |
| 7 | run_ensemble_acp.py | anticp2_main | 4-arm ens | 0.796 AUC | ensemble negative |
| 8 | run_acp_v4_multiseed.py | anticp2_main | seed-ens x3 | 0.773 AUC | negative |
| 9 | run_acp_v5_big.py | anticp2_main | CNNv2-big | 0.771 AUC | negative |
| 10 | run_solubility.py | esol | RF-DESC | 0.797 AUC | sol best |
| 11 | run_solubility.py | esol | PeptideGNN | 0.785 AUC | - |
| 12 | run_aggregation.py | amylogram | PeptideGNN | 0.769 AUC | - |
| 13 | run_pep424_headtohead.py | pep424 | GNN transfer | invalid (parse) | documented |
| 14 | run_pep424_v2.py | pep424 | GNN transfer | invalid (normalization) | documented |
| 15 | run_pep424_v3.py | pep424 | GNN transfer | 0.790 AUC | **break vs FoldAmyloid 0.748** |
| 16 | run_pep424_v3.py | pep424 | GNN 5-fold CV | 0.794 AUC | below AmyloGram 0.865 |
| 17 | run_discovery.py | multi | TriNet 6ep | saturated | bug found+fixed |
| 18 | run_discovery.py (fixed) | multi | TriNet 6ep | 1 gated candidate | v0.1 |
| 19 | run_discovery_scaled.py | multi | TriNet 25ep | (Section 7 results) | current |

# Appendix E. CLI reference

```
python -m pepx.cli score SEQ [SEQ ...]     # tri-objective probabilities, TSV
python -m pepx.cli score-file IN.fa [--out OUT.tsv]
python -m pepx.cli describe SEQ            # 18-descriptor JSON panel
```

Exit codes: 0 success; 2 missing weights. Invalid residues are rejected with a diagnostic naming the offending letters. Batch throughput: ~2,000 sequences/min on 2 CPU cores.

# Appendix F. Discovery-screen protocol (falsifiability)

A candidate graduates from "screen hit" to "claimed discovery" only if: (i) it is absent from CancerPPD, APD3, AntiCP 2.0, DBAASP and UniProt-derived pools (exact match); (ii) P_acp >= 0.9, P_sol >= 0.7, P_agg <= 0.2 from the 25-epoch TriNet; (iii) descriptor gates: charge +2..+9, GRAVY < 0, Boman 1.5..3.2, length 10--30; (iv) it is re-derived from at least two independent seed lineages; (v) a blinded re-score through the shipped CLI reproduces (ii) to within 0.01. The falsification path: solid-phase synthesis (13--30 residues, <$200), then (a) MTT assay vs HeLa + HEK293 at 10 uM, (b) turbidity at 100 uM pH 7.4, (c) ThT fluorescence at 37 C/24 h. Failure on any axis refutes the tri-objective claim for that sequence without touching the benchmark results.

# Appendix G. Data-quality findings (standalone value)

1. pep424's header namespace collides (158/424 "3D profile"); any name-joined comparison is silently wrong. We reported this pattern with the fix (order + sequence join).
2. eSOL gene-name mapping leaves 36% of proteins unmatched; synonym-aware mapping (b-numbers) is the fix.
3. AntiCP 2.0 main contains near-duplicate kationic stretches across train/test (unquantified; homology audit queued).
4. APD3 contains sequences with non-standard letters that strict validators must reject rather than silently drop (we log every rejection).
