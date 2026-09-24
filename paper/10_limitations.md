# Limitations

We list every limitation we are aware of, in decreasing order of consequence. A result section that does not name its limits is marketing; this project exists to be believed.

## Conformational blindness

pepx scores sequences, not structures. The PDB analysis (Section 17) shows the canonical ACPs are micelle-induced alpha-helices; nothing in our models represents induced fit, membrane partitioning depth, or orientation. The Chou-Fasman scales (Section 2) are a 1978 statistical proxy. A candidate that passes all gates can still fail because it does not fold correctly in a membrane-mimetic environment. This is the strongest argument for the falsification protocol (Appendix F) being wet-lab-first rather than model-forever.

## Label provenance and noise

AntiCP 2.0's negatives are "non-ACP AMPs" (main) and random UniProt-derived peptides (alternate) - negative by absence of annotation, not by experimental falsification. Some fraction are undiscovered ACPs; that fraction is a hard ceiling on achievable AUROC, and our convergence of ten-plus configurations at 0.77-0.80 against the published 0.83 is consistent with hitting it (Section 5). We state the 0.83 bar as *published but not independently reproduced by us*; the AntiCP 2.0 paper's own protocol details (exact negative pool, exact preprocessing) are not fully recoverable from the manuscript, which is itself a reproducibility finding we record here.

## Aggregation ground truth is assay-bound

WALTZ-DB/pep424 labels come from specific assays (ThT, EM, FTIR) at specific concentrations and pH. A sequence labeled non-amyloid at 100 uM pH 7 may aggregate at 500 uM pH 5. Our 0.8391 transfer result is therefore a statement about *this label function*, not about aggregation as a physical property. The AmyloGram 0.865 CV bar standing unbeaten (0.7947) is consistent with within-assay learnability being higher than cross-assay transfer.

## Single-seed headline numbers

Headline AUROCs are single locked-seed evaluations. Where we measured spread (10-config ACP sweep 0.771-0.801; 12-trial HPO test-at-best-val 0.731-0.788), within-task variance is 0.02-0.05 AUROC. The pep424 break (+9.1) is far outside this band; the AntiCP2 gap (-2.7) is inside it. Confidence language in this paper follows that arithmetic: the break is robust, the gap is directional.

## Discovery candidates are computational

FEKEAKKIEIKRH has never touched a cell. Every claim about it is a model claim plus a descriptor-gate claim plus a novelty screen against 203 external sequences (0 hits). The anionic drift of the 25-epoch screen (Section 15) demonstrates that these models can be confidently wrong in regions far from training support. We make no therapeutic claim.

## Scale

2 CPU cores, 1.9 GB RAM, no GPU. Larger models, longer searches, k-fold ensembling at scale, and learned embeddings (ESM-class) were out of reach and are named as future work, not attempted and quietly dropped.
