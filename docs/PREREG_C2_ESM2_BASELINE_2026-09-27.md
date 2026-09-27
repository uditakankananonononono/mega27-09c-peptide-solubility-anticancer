# PREREGISTRATION C2 - ESM-2 embedding baseline (locked 2026-09-27, BEFORE
# any C2 embedding or scoring)

Verdict item C2: "ESM-2/ProtT5 embedding baseline on existing splits."
Scope locked to ESM-2 (ProtT5 does not fit this 2GB sandbox; disclosed).

## Model (locked; same recipe as 09b B1 arm2, cross-repo consistency disclosed)
- Frozen esm2_t12_35M_UR50D, final-layer mean-pooled 480-d embeddings.
- Head: logistic regression (lbfgs, C=1.0, max_iter=2000) per task, trained
  on that task's training split ONLY. Single scoring per task on its
  evaluation split. No threshold tuning (0.5 for MCC).

## Splits and comparators (all pre-existing; nothing re-split)
1. anticp2_main (PRIMARY): official train/validation partition (verified
   line-identical to the paper's files, COMPARATOR_ANTICP2.md). Comparators,
   verbatim published validation rows: AUROC 0.83, MCC 0.51. Decision:
   BEAT if metric > comparator AND 10k bootstrap lower95 > comparator;
   MATCH if CI spans; MISS otherwise. Both metrics reported as-is.
2. anticp2_alternate: official partition; REPORT-ONLY (no verbatim published
   alternate validation rows in the repo - disclosed).
3. esol seed-7 split (results/esol_split_seed7.json committed indices,
   30% solubility cutoff): REPORT-ONLY baseline (no external comparator on
   this exact split).
4. pep424 aggregation: train on amylogram_full hexapeptides; score pep424
   by hexapeptide-window max-pool (the repo's existing convention,
   run_pep424_headtohead.py), canonical-scorable subset (n=261, 163
   non-canonical excluded - disclosed). Comparators: AmyloGram published
   AUC 0.865 on the FULL 424 (different n - disclosed caveat), and
   FoldAmyloid same-subset AUC 0.5688 (committed pep424_headtohead.json) as
   the same-partition row. BEAT/MATCH/MISS vs the same-partition row only.

## Falsifier (locked)
- Label-shuffle head on anticp2_main (seed-13 permutation of train labels):
  validation AUROC must fall in [0.40, 0.60]; otherwise void until fixed.

## Reporting
results/c2_esm2_baseline.json, all four tasks, as-is.
