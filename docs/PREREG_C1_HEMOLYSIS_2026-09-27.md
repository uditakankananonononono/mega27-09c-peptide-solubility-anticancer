# PREREGISTRATION C1 - Hemolysis fifth task (locked 2026-09-27, BEFORE any
# hemolysis training or scoring)

Verdict item C1: "Hemolysis classifier (HemoPI-class dataset) as a fifth
task."

## Data (source of record)
- HemoPI official dataset files (Chaudhary et al., Sci Rep 2016, srep22843,
  PMC4782144), retrieved 2026-09-27 from
  https://webs.iiitd.edu.in/raghava/hemopi/data/HemoPI_{1,2,3}_dataset/{main,validation}/{pos,neg}.fa
  sha256 per file in data/raw/hemopi/SHA256SUMS.txt (committed).
- Splits: the published main/validation partitions, used exactly as shipped.
  HemoPI-1: 442+442 main, 110+110 validation. HemoPI-2: 442+370 main,
  110+92 validation. HemoPI-3: 708+590 main, 177+148 validation.

## Comparators (published validation rows, Table 5 of PMC4782144, verbatim)
- HemoPI-1 validation (hybrid motif+SVM, the paper's best): Sn 96.4,
  Sp 99.1, Acc 96.4, MCC 0.93.
- HemoPI-2 validation: Sn 78.2, Sp 78.3, Acc 75.7, MCC 0.51.

## Model and training (locked)
- PeptideCNNv2 (the repo's shipped architecture) + train_single_task with
  the exact AntiCP2 recipe: epochs 40, batch 64, lr 1e-3, weight decay 1e-4,
  val_frac 0.15 (class-stratified, seed 7), patience 8, wide=True,
  max_len 60, seed 7.
- NO ACP-specific augmentation for hemolysis (augmented_acp_train is
  ACP-only) - disclosed difference from the AntiCP2-main recipe.
- Trained on main split ONLY; scored ONCE on the validation split per
  dataset. No threshold tuning on validation (threshold 0.5, paper
  convention).

## Endpoints and decision rules (locked)
- Primary endpoint: HemoPI-1 validation MCC.
  BEAT: MCC > 0.93 AND 10,000-replicate bootstrap lower95 of our MCC > 0.93.
  MATCH: bootstrap CI spans 0.93. MISS: upper95 < 0.93. All reported as-is.
- Secondary endpoint: HemoPI-2 validation MCC vs 0.51 (same rule).
- Tertiary, report-only: HemoPI-1/2 validation AUROC; HemoPI-3 validation
  MCC + AUROC (paper reports HemoPI-3 CV only - no published validation row,
  so no comparator gate; disclosed).

## Falsifier (locked)
- A label-shuffle model (seed-13 permutation of main-split labels, identical
  architecture/recipe) must score |MCC| < 0.15 on HemoPI-1 validation; if
  not, the harness is void until fixed.

## Reporting
results/c1_hemopi.json, per-dataset metrics + CIs + comparator deltas,
as-is.
