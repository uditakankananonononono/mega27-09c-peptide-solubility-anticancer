# Architecture and training protocol, in full

Everything needed to reimplement each model from this document alone. Tensors flow batch-first; all models end in a single logit per task.

## PeptideCNN (baseline deep model)

Embedding $E \in \mathbb{R}^{20 \times 32}$ (padding index 0 shares the table with no residue; sequences right-padded to 60). Three parallel Conv1d branches, kernels 3/5/7, 64 filters each, ReLU, global max-pool per branch. Concatenated $3 \times 64 = 192$ features with the 18-descriptor vector (train-statistics standardized), then Dropout(0.3) - Linear(210, 128) - ReLU - Linear(128, 1). Parameter count: 45,153 (computed, `sum(p.numel())`). This is the "cheap deep baseline" every classical model must beat to justify depth; on AntiCP 2.0 main it does not (0.78-0.79 vs ET-DPC 0.799), which is part of why we trust the classical numbers.

## PeptideGNN

Nodes are residues; edges connect sequence-adjacent pairs (both directions, no self-loops). Node features: per-residue one-hot (20) concatenated with per-residue physicochemical scales (hydrophobicity, charge contribution, mass, Chou-Fasman pair) - the mass feature is why the P/Q/V/W/Y table scramble moved the pep424 number by 4.9 AUROC points, and why the fix was load-bearing. Two graph convolution layers (mean aggregation, 48 hidden), global mean+max pool concatenated, MLP head. Parameter count: 16,481 (computed). No edge features, no attention: the chain graph is the inductive bias, and on 419-sequence pep424 it is enough to beat a packing-density published tool by +9.1.

## PeptideCNNv2 (AntiCP2 challenger)

Residual dilated convolutions (dilation 1, 2, 4; 96 channels; kernel 3; pre-activation residual blocks), additive attention pooling (Section 19.3 with the padding-mask caveat that cost us early runs), fused with a 438-wide feature block (18 descriptors + 20 AAC + 400 DPC), head Dropout(0.35) - Linear(2*96+438, 128) - ReLU - Linear(128,1). Parameter count: 253,746 (computed). Best main-split single model: 0.8011; HPO-chosen variant: 0.8029 (Section 5.4).

## TriNet (the discovery machine)

Shared trunk: embedding 48, parallel Conv1d 3/5/7 (96 filters each), global max-pool, fused with the 18-descriptor block, Dropout(0.3) - Linear(306, 128) - ReLU. Three heads, each Dropout(0.3) - Linear(128, 1): acp / sol / agg. Parameter count: 110,099 (computed). Trained round-robin over the three task loaders, one epoch each per cycle, class-weighted BCE (pos_weight = neg/pos per task batch), Adam 1e-3 with weight decay 1e-4, 6 cycles (v0.1, shipped) and 25 (the drifted negative). Held-out AUROCs at ship time: acp 0.79 / sol 0.80 / agg 0.78 - deliberately mid-strength heads: the screen relies on gates + ranking, not on any head being a frontier classifier.

## Training protocol common to all runs

Seed 7 everywhere (torch, numpy, python). Stratified validation carve 15% of train (class-balanced indices), early stopping on validation AUROC, patience 8 (5 for HPO trials), batch 64, max lengths 60 (peptide tasks) / 300 (eSOL full proteins) / 32 (amyloid hexapeptides). Descriptor standardization: train-split statistics only, per task - the single most important implementation rule in the repo; violating it once saturated every score to 0/1 (documented in the git history of the discovery screen). Loss: BCEWithLogits with batch-computed pos_weight for imbalanced tasks (eSOL 1.65:1, amyloid 2.5:1).

## Evaluation protocol

Locked test sets are touched exactly once per run, after early stopping selects the checkpoint. AUROC is the headline metric everywhere (threshold-free, matches published bars); accuracy/MCC at 0.5 accompany it. The 10-config sweep and the 12-trial HPO record test-at-selection for every configuration - the spread is published in this paper (0.731-0.788) precisely so that no one, including us, can mistake a lucky trial for progress.
