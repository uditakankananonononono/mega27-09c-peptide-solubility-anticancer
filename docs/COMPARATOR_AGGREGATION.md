# Comparator documentation - aggregation arm (locked protocol C1, task C)

## Sources (pinned)
- PASTA 2.0 paper full text: Walsh et al., NAR 2014 (gku399), PMC4086119,
  data/external/pasta2_pmc4086119.xml (NCBI efetch, sha256 recorded in git).
- AmyloGram pep424 benchmark: michbur/AmyloGramAnalysis (already pinned in
  data/raw/aggregation/, incl. FoldAmyloid same-partition predictions).

## Verified published rows (PASTA 2.0 paper, verbatim)
- Total AUC 85.73 (= 0.8573) on their validation benchmark; FoldAmyloid next
  best, AUC 2.42 points worse (= 0.8331).
- 90%-specificity threshold: Sen 30.24, Spc 90.00, Acc 80.23, MCC 0.22.
- 85%-specificity threshold: Sen 40.87, Spc 84.95, Acc 77.77, MCC 0.24.

## Partition caveat (honest)
PASTA 2.0's published AUC is on THEIR validation benchmark, not on pep424.
Identical-partition PASTA 2.0 scores require running the tool on the 419
aligned pep424 sequences: standalone is behind a registration form
(biocomputingup.it) and the webserver needs a browser pass - QUEUED for the
next browser token, not fabricated.

## Same-partition state NOW (pep424, n=419 aligned)
- FoldAmyloid (same-partition predictions from the AmyloGram analysis repo):
  AUC 0.7480 (results/pep424_v3.json).
- Our transfer pepx-GNN (trained on AmyloGram hexapeptides, scored on pep424):
  AUC 0.8391 - BEATS FoldAmyloid on the identical partition (+0.091).
- Our 5-fold CV on pep424: AUC 0.7947.
- AmyloGram published CV AUC 0.865 (cross-partition number - their own CV).
Verdict: one same-partition comparator beaten (FoldAmyloid); strongest
comparators (PASTA 2.0, AmyloGram) not yet evaluated on the identical
partition - PASTA 2.0 webserver run queued; per C1 the arm's beat gate
requires the strongest comparator on the identical partition with bootstrap CI.
