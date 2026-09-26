# Pre-registration addendum - 2026-09-27 (09c ACP arm, improvement iteration 1)
Locked BEFORE any stacked-model scoring. Prior registrations unchanged.

## Target
Close the honest gap documented in docs/COMPARATOR_ANTICP2.md: published
AntiCP 2.0 ETree-DPC validation AUROC 0.83 / MCC 0.51 on the OFFICIAL main
validation set (verified identical partition) vs our cnnv2-hpo-best
AUROC 0.80286 / MCC 0.4638.

## Locked iteration-1 protocol
Stacked generalization, train-side only:
- Base learners trained on the official train split with out-of-fold (5-fold,
  stratified, seed 2709) probability outputs:
  (a) CNNv2 at the committed HPO config (emb 32, ch 96, blocks 2,
      dropout 0.2781, lr 0.001729);
  (b) ET-DPC replication (ExtraTrees 400 trees on dipeptide composition);
  (c) AAC + train-mined motif-presence logistic (MERCI-style: motifs mined
      from the train fold only, never from test).
- Meta-learner: logistic regression on the 3 out-of-fold probability columns.
- ONE evaluation on the official validation set (n=344). No test-informed
  selection; if the stack loses to cnnv2-hpo-best alone on validation AUC,
  that is reported honestly as a negative iteration.
- Comparison vs published: 10,000-replicate bootstrap on the official
  validation predictions; BEAT = stacked AUROC > 0.83 AND MCC > 0.51 with
  the AUROC bootstrap CI of (stack - published point value) excluding 0
  (published rows are point values from the paper; their variance is not
  recoverable from the paper - disclosed limitation).
- All artifacts (OOF predictions, test predictions, bootstrap output)
  committed under results/ with this addendum referenced.
