# Locked comparator updates: AntiCP and DeepSol

## AntiCP 2.0 official-partition identity

The science branch pinned an official AntiCP 2.0 dataset and paper in `docs/COMPARATOR_ANTICP2.md`. That document reports a split-identity check (reported as a zero-line diff) against the repository split. The published benchmark values used in the locked gate are AUROC 0.83 and MCC 0.51 for its named comparator row. Our prior CNNv2 HPO-best row (reported in the comparator document) is AUROC 0.8029 and MCC 0.4638 on that official test. We do not relabel these as a win. A published row is still a fixed point-value comparison, not a paired same-code rerun unless predictions are aligned and available.

## Two attempts, two honest negatives

The first iteration was locked before scoring by the addendum in the science branch. The first committed stack result records an OOF meta-logistic stack of CNNv2-HPO, ET-DPC and motif logistic components. Training-internal OOF AUROC was 0.8532, but the single official test achieved AUROC **0.7884**, MCC **0.4192**, below both the published 0.83/0.51 gate and our earlier CNNv2 benchmark row. The bootstrap difference against the *published point AUROC* had an interval spanning zero. This is a model-selection/transfer failure, not a near-win to round up.

Iteration 2 used an a priori equal-weight strategy selected using train-internal CV before the one official-test score. The second committed stack result records the selected seven-seed CNN plus ET-DPC combination, test AUROC **0.8121** and MCC **0.4887**. It improved on the earlier internal CNNv2-HPO row but remains below the two-part published gate; `beat` is false. Its bootstrap difference against the published point AUROC also spans zero. It would be incorrect to call an interval crossing zero proof of equivalence, or to search unregistered extra weights on the same test set until a favorable number appears. The source addendum says to consult a redirection path after this second negative.

| Attempt | AUROC | MCC | Gate |
|:--|--:|--:|:--|
| Prior CNNv2-HPO | 0.8029 | 0.4638 | No |
| OOF stack, iteration 1 | 0.7884 | 0.4192 | No |
| CNN + ET-DPC, iteration 2 | 0.8121 | 0.4887 | No |
| Published comparator | 0.8300 | 0.5100 | Reference |

The OOF stack had train-internal AUROC 0.8532. The second combination was selected on train-internal CV. Neither internal selection figure belongs in the official-test columns.

The table intentionally separates train-internal selection and official-test evaluation. It does not promote repeated exposure to the same test partition as an independent validation exercise.

## DeepSol and PASTA 2.0: protocols versus results

`docs/COMPARATOR_AGGREGATION.md` pins a PASTA 2.0 paper comparator with a stated AUROC 0.8573 and MCC 0.22-0.24, while noting a partition mismatch. The project's 0.8391 versus FoldAmyloid 0.7480 is an aligned pep424 comparison, but there is no committed same-partition PASTA 2.0 run to add to that win. [PENDING: a lawful, access-compatible same-partition PASTA 2.0 evaluation; do not treat literature point values as a paired score.] The locked DeepSol arm has a pinned official split and predeclared models according to its addendum; no final scored result has been committed in this paper update. [PENDING: the official DeepSol final test and bootstrap result.] Intermediate epochs or validation numbers are not submission-ready evidence of a test win.

## Consequences for claims

AntiCP remains an honest negative, not a hidden step along a guaranteed progress curve. The pep424 aligned FoldAmyloid result remains separate, and the multi-objective screen remains computational. None of the AntiCP iteration numbers proves differential cancer-cell selectivity or clinical benefit; those require independent assays. The result is the documented limit of these two locked improvement approaches on one official partition, with clear next-work conditions and without changing benchmark definitions after inspection.
