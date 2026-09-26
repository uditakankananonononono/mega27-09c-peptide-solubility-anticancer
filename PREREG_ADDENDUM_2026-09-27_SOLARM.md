# Pre-registration addendum - 2026-09-27 (09c SOLUBILITY arm, DeepSol same-split)
Locked BEFORE any model is scored on the DeepSol partition. Prior registrations unchanged.

## Data (pinned, sha256 in data/raw/deepsol/SHA256SUMS.txt)
DeepSol official deposit, Zenodo DOI 10.5281/zenodo.1162886
(sameerkhurana10/DSOL_rv0.2-v0.3.zip), files data/{train,val,test}_{src,tgt}:
protein sequences + binary labels on DeepSol's OWN partition -
train 62,478 / val 6,942 / independent test 1,999. This is the identical
partition the published numbers come from.

## Published comparator numbers (verified from the paper record)
Khurana et al., Bioinformatics 34(15):2605-2613 (2018), PMC6355112,
DOI 10.1093/bioinformatics/bty166. Abstract (verbatim): DeepSol "attained an
accuracy of 0.77 and Matthew's correlation coefficient of 0.55" on the
independent test set.

## Locked protocol
- Models: (a) ET-DPC (ExtraTrees 400, dipeptide composition); (b) CNNv2 at the
  committed HPO config adapted for longer proteins (max_len 512);
  (c) a priori equal-weight mean of (a)+(b) probabilities (no learned combiner -
  ACP-arm lesson).
- Model selection/early stopping on DeepSol's OWN val split only. ONE
  evaluation per model on the independent test (1,999). No test-informed
  iteration; any later change requires a new appended addendum.
- Metrics: accuracy, MCC (published primary pair), plus AUC reported
  secondarily. BEAT = accuracy > 0.77 AND MCC > 0.55 on the identical test
  partition with 10,000-replicate bootstrap CI of (ours - published) excluding
  0 for both metrics (published rows are point values - disclosed limitation).
- Honest-negative clause: negative = reported as-is and arm goes to rule-6
  redirection consult.
