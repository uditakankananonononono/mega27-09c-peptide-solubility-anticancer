# 09c research/data tool audit - 2026-09-25

Strict candidate count: 24, **not 40**. A candidate needs execution evidence, not merely a citation. Results were committed before this audit and have not all been rerun here. Count one package/tool once even if used repeatedly. The table is a traceable candidate ledger, not a declaration that all 24 were independently reproduced.

| Candidate | Research use | Code / committed output |
|---|---|---|
| PyTorch | CNN/GNN/TriNet fitting | `src/pepx/models.py`, `results/trinet.pt` |
| scikit-learn | classifiers, splits, metrics | `runs/tool_analyses.py`, `results/baselines_anticp2.json` |
| NumPy | descriptor arrays, numerical work | `src/pepx/encoders.py`, `results/descriptor_crossvalidation.json` |
| SciPy | statistics and comparisons | `runs/tool_analyses.py`, `results/tool_analyses.json` |
| pandas | benchmark tables | `src/pepx/datasets.py`, `results/dataset_stats.json` |
| modlAMP | independent descriptor check | `results/descriptor_crossvalidation.json` (historical output; code invocation not retained, provisional) |
| peptides.py | independent descriptor check | `results/descriptor_crossvalidation.json` (historical output; code invocation not retained, provisional) |
| matplotlib | benchmark plots | `runs/tool_figures.py`, `paper/fig5_umap.png` |
| XGBoost | ACP baseline | `runs/tool_analyses.py`, `results/tool_analyses.json` |
| LightGBM | ACP baseline | `runs/tool_analyses.py`, `results/tool_analyses.json` |
| statsmodels | logistic baseline | `runs/tool_analyses.py`, `results/tool_analyses.json` |
| NetworkX | similarity graph | `runs/tool_analyses.py`, `results/tool_analyses.json` |
| propy3 | independent CTD descriptors | `runs/tool_analyses.py`, `results/tool_analyses.json` |
| pyteomics | independent peptide mass check | `runs/tool_analyses.py`, `results/tool_analyses.json` |
| biotite | local alignment and PDB metrics | `runs/tool_analyses.py`, `results/tool_analyses.json` |
| UMAP | descriptor embedding | `runs/tool_figures.py`, `results/tool_figures.json` |
| SHAP | solubility attribution | `runs/tool_figures.py`, `results/tool_figures.json` |
| logomaker | ACP logos | `runs/tool_figures.py`, `results/tool_figures.json` |
| upsetplot | overlap visualization | `runs/tool_figures.py`, `results/tool_figures.json` |
| seaborn | statistical plots | `runs/tool_figures.py`, `paper/fig5_umap.png` |
| numba | JIT kernel experiment | `runs/tool_engineering.py`, `results/tool_engineering.json` |
| onnxruntime | independent inference parity | `runs/tool_engineering.py`, `results/tool_engineering.json` |
| imbalanced-learn | SMOTE comparison | `runs/tool_engineering.py`, `results/tool_engineering.json` |
| Optuna | CNNv2 HPO | `runs/run_optuna_acp.py`, `results/optuna_cnnv2_anticp2.json` |

Excluded: Biopython's static DIWV table (data source, not executed research library); AntiCP 2.0, CancerPPD, APD3, eSOL, UniProt, AmyloGram, DBAASP, RCSB PDB, ChEMBL (data/literature sources); FoldAmyloid prediction table (comparator, not tool executed); ruff, mypy, pytest-cov, hypothesis (QA/build tooling); ONNX export format (not a separate tool from onnxruntime); FastAPI/pydantic (serving infrastructure, not research); pytest/git/pandoc/curl. If any candidate lacks replayable evidence, subtract it until it can be confirmed. The 40 gate is open.
