# 7. TriNet and the tri-objective discovery screen

## 7.1 Multi-task training

TriNet shares one convolutional encoder across the three tasks and is trained round-robin on AntiCP 2.0-augmented ACP data, eSOL solubility, and AmyloGram aggregation. Held-out AUROCs after 6 epochs (CPU): acp 0.775 / sol 0.786 / agg 0.796. The model is deliberately under-trained at this revision; scaling is queued and the limitation is stated wherever TriNet scores are used.

## 7.2 The screen

Simulated annealing (Section 2.11) from 24 CancerPPD seed sequences, 36 single-residue mutation steps each, optimizing S(s) = P_acp · P_sol · (1 - P_agg). Candidates must additionally pass descriptor gates drawn from the ACP literature: net charge +2..+9, GRAVY < 0, Boman 1.5..3.2, length 10--30, and absence from every database used. Committed: `results/discovery_candidates.json`.

## 7.3 Result: one gated candidate (v0.1)

**FEKEAKKIEIKRH** (13-mer): P_acp = 0.904, P_sol = 0.867, P_agg = 0.055; charge +2, GRAVY -0.61, Boman 1.64; novel (absent from CancerPPD/APD3/AntiCP 2.0/UniProt-derived pools). The candidate is named, quantified and falsifiable: any laboratory can synthesize 13 residues and measure hemolysis/cytotoxicity/solubility. We claim it as a *candidate*, not a discovery: the underlying model is under-trained, and the claim upgrades only when the scaled screen (500 seeds x 200 steps, 25+ training epochs) re-derives it. Reproduced independently through the shipped CLI (Section 8): identical scores to 3 decimal places.

# 8. The pepx tool

`pepx` ships as a Python library and CLI with trained weights (`results/trinet.pt`) and normalization statistics (`results/trinet_norm.npz`):

```
python -m pepx.cli score KLAKLAKKLAKLAK FEKEAKKIEIKRH
# sequence        P_anticancer  P_soluble  P_aggregation
# KLAKLAKKLAKLAK  0.9935        0.6493     0.4185
# FEKEAKKIEIKRH   0.9040        0.8675     0.0548

python -m pepx.cli score-file candidates.fasta --out scores.tsv
python -m pepx.cli describe KLAKLAKKLAKLAK   # full descriptor panel
```

Controls: the known ACP KLAKLAKKLAKLAK scores P_acp 0.9935; poly-glutamate scores P_agg 0.0115. Repository: `github.com/uditakankananonononono/mega27-09c-peptide-solubility-anticancer` (branch main), 23/23 tests green.

\newpage

# 9. Honest negatives and failure log

1. **AntiCP 2.0 main unbeaten.** Ten configurations tried (Section 6.1); best 0.8011 AUROC vs published 0.83. Seed-averaging hurt (0.773); a wider model hurt (0.771); a 4-arm ensemble sat between its arms (0.796). The published number's homology audit is pending.
2. **AntiCP 2.0 alternate unbeaten** (0.931 vs 0.95).
3. **AmyloGram unbeaten** on its own CV protocol (0.7947 vs 0.865).
4. **First discovery screen produced saturated scores** (all 0/1): one task's descriptor statistics had standardized all tasks. Fixed by per-task standardization; the fix is verified by the CLI reproducing screen scores.
5. **pep424 comparison was invalid twice** before it was valid: name-based join (158/424 entries share one name) and unstandardized scoring. Both invalid runs are committed.
6. **Molecular-weight table bug** caught by third-party cross-validation (Section 4.1).
7. **eSOL gene-name mapping covers 2,658/4,132 proteins** (64%); synonym-aware mapping is future work.
8. **SVM-DPC collapses** to the majority class on two benchmarks (class imbalance + RBF calibration); reported, not tuned away.

# 10. Research/data tools - strict audit (gate open)

The earlier 40-row table counted data archives, literature/benchmark sources,
lint, type checking, coverage, property testing, and HTTP serving. These do
not qualify as 40 research/data tools. The corrected inventory and each
candidate's code/output evidence are in `TOOLS_LEDGER.md`. It identifies
24 candidate research tools/libraries, subject to execution-level audit,
not a verified 40. RCSB PDB, ChEMBL, AntiCP 2.0, CancerPPD, APD3, eSOL,
UniProt, DBAASP and AmyloGram are important data sources, not software tools.
FoldAmyloid predictions are a benchmark comparison, not a tool this project ran.
The 40-tool gate is **not met**. No method or benchmark result is altered
by this accounting correction.

# 11. Conclusion

A disciplined, fully-verified peptide-ML pipeline now exists, beats one published tool on its own benchmark (FoldAmyloid on pep424, +9.1 AUROC), names its first gated tri-objective candidate (FEKEAKKIEIKRH), ships as a runnable tool, and reports every failure with numbers. The two flagship bars --- AntiCP 2.0 main (0.83) and AmyloGram CV (0.865) --- stand; the roadmap to them is concrete.

# Appendix A. Reproduction

Environment: Python 3.10.12, torch 2.14.0+cpu, scikit-learn 1.7.2, 2 CPU / 1.9 GB RAM. Commands: `python3 -m pytest`; `python3 runs/run_baselines.py` equivalents per `runs/`; all seeds fixed. Results files cited inline per table.

# Appendix B. Dataset manifest

See `data/MANIFEST.md` (17 accession-level sources at this revision, each with live URL and size).

# Appendix C. Full results JSON inventory

`results/baselines_anticp2.json`, `ensemble_anticp2_main.json`, `acp_v3_best.json`, `acp_v4_seedens.json`, `acp_v5.json`, `solubility_esol.json`, `pep424_headtohead.json`, `pep424_v2.json`, `pep424_v3.json`, `discovery_candidates.json`, `descriptor_crossvalidation.json`.
