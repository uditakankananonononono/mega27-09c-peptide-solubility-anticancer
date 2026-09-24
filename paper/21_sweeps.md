# The AntiCP 2.0 main-split campaign, unabridged

The 0.83 bar was attacked ten-plus times; the complete committed trail (all numbers from the cited results files, locked test, seed 7):

| attempt | model / idea | test AUROC | lesson |
|---|---|---|---|
| baselines | SVM-AAC | 0.579 | composition alone is weak |
| baselines | RF-DPC | 0.779 | motifs carry the signal |
| baselines | ET-DPC | 0.799 | best classical |
| v2 | CNNv2 + APD3 augmentation | 0.781 | augmentation neutral (`cnnv2_aug_anticp2_main.json`) |
| v3 | CNNv2 dropout/lr variants | 0.788 | tuning moves +-0.01 (`acp_v3_best.json`) |
| ensemble | ET+CNNv2+GNN+CNN | 0.795 | ensemble of uncorrelated models *underperforms* its best member (`ensemble_anticp2_main.json`) |
| v4 | CNNv2 seed-ensemble (7/17/27) | 0.773 | seed-averaging hurt (`acp_v4_seedens.json`) |
| v5 | CNNv2-big + ET blend | 0.771 | capacity and blending both negative (`acp_v5.json`) |
| v6 | CNNv2 (40-epoch, tuned) | 0.801 | best single deep model |
| HPO | Optuna 12-trial TPE | 0.803 | search converges to the same ceiling (`optuna_cnnv2_anticp2.json`) |

Three post-mortems worth more than the numbers. **(i) The ensemble negative.** Averaging four models with genuinely different inductive biases *lost* 0.6 points to the best member - the errors were not independent: every model had learned the same lysine-richness direction (Figure 7 makes the shared signal visible), so averaging averaged away what little orthogonal signal existed. **(ii) The seed-ensemble negative.** Three seeds, same architecture, averaged probabilities: worse than the median single seed. With small data the seed variance is in *which* near-homologs get memorized (Section 18 audit), and averaging correlates the memorization rather than canceling it. **(iii) The HPO plateau.** Bayesian search over five dimensions landed 0.2 points from hand tuning. When three very different optimization pressures (manual, ensembling, Bayesian) converge on 0.77-0.80, the bound is in the data, and the Section 18 leakage audit says precisely where: generalization beyond the near-homolog stratum is the missing piece, and no amount of model-side effort recovers information the labels do not contain.

## The alternate split, briefly

Every classical model scores higher on the alternate split (ET-DPC 0.931 vs published 0.95) than on main - because its negatives are random UniProt fragments rather than non-ACP AMPs. The 1.9-point gap to the published bar is the same story at smaller magnitude; we spent our effort on the main split because its negative pool makes it the harder and more honest task.
