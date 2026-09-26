# Comparator documentation - AntiCP 2.0 (locked protocol C1, task A)

## Sources (pinned, sha256 in data/external/anticp2_official/SHA256SUMS.txt)
- Paper PDF: Agrawal et al., "AntiCP 2.0: an updated model for predicting
  anticancer peptides", Brief Bioinform 2021 (bbaa153), fetched from the
  authors' lab reprint https://webs.iiitd.edu.in/raghava/reprints/AntiCP2.pdf
- Official dataset files: https://webs.iiitd.edu.in/raghava/anticp2/download.php
  (add_data/{pos,neg}_{train,test}_{main,alternate})

## Split-identity verification (LIVE, 2026-09-27)
Sorted-content diff of repo data/raw/anticp2/{pos,neg}_{train,test}_main vs the
official download = 0 differing lines for all four files. Our anticp2_main split
IS the paper's official train/validation partition. Same-split comparison is
therefore valid without re-splitting.

## Published rows (verbatim from paper PDF)
Main dataset, best model = DPC-based ExtraTrees (Ntree=400), Table 2:
- Training (5-fold CV): Sen 74.06, Spc 76.52, Acc 75.29, MCC 0.51, AUROC 0.83
- Validation (official test): Sen 77.46, Spc 73.41, Acc 75.43, MCC 0.51, AUROC 0.83
Table 9 (existing-method comparison on the main validation set):
- AntiCP 2.0 MCC 0.47; ACPred-Fuse MCC 0.32 (text: "AntiCP 2.0 got MCC 0.47,
  whereas ACPred-Fuse got MCC 0.32")
Abstract best: "MCC of 0.51 and 0.83 AUROC on the training dataset" (main).

## Our committed numbers on the IDENTICAL official validation set (n=344)
- PeptideCNNv2+aug (results/cnnv2_aug_anticp2_main.json, commit history):
  AUC 0.8029, MCC 0.4638
- ENSEMBLE(ET+CNNv2+GNN+CNN) (results/ensemble_anticp2_main.json):
  AUC 0.7958, MCC 0.4347
- ET-DPC in-house replication of the paper's best model
  (results/baselines_anticp2.json): AUC 0.7990, MCC 0.4467

## Honest verdict (2026-09-27): NOT BEATING
Primary metric gap on the identical partition: AUC 0.8029 vs 0.83 (-0.027),
MCC 0.4638 vs 0.51 (-0.046). Our ET-DPC replication (0.799/0.4467) also lands
below the paper's published validation rows, consistent with implementation
differences; the comparison target remains the published numbers per C1.
Next: improvement iterations, each locked as a new prereg addendum before
scoring (feature/architecture/ensemble-direction), until the primary metric
beats 0.83 AUROC / 0.51 MCC on this exact partition with a 10,000-replicate
bootstrap CI excluding 0, or a rule-6 pivot redirects the arm.
