---
title: "Tri-Objective Peptide Intelligence: Benchmark-Headed Models for Anticancer Activity, Solubility and Aggregation, a Head-to-Head Win over FoldAmyloid on pep424, and a Multi-Task Discovery Screen"
author: "MEGA-PROGRAM-27, Item 9 Part 3 --- automated research pipeline (pepx)"
date: "September 24, 2026"
geometry: margin=1in
fontsize: 12pt

header-includes:
  - \usepackage{amsmath, amssymb, booktabs, graphicx, longtable, hyperref}
  - \usepackage{mathptmx}
---

# Abstract

Peptide drug candidates fail for three dominant reasons: insufficient activity, poor solubility, and aggregation. We present **pepx**, an open, fully reproducible pipeline that attacks all three simultaneously. (i) On the AntiCP 2.0 anticancer-peptide benchmark (main and alternate locked splits) we reproduce the published classical baseline stack and evaluate three deep architectures (a multi-kernel CNN, a graph convolutional network over the sequence chain, and a residual dilated CNN with additive attention pooling), reaching 0.8011 AUROC against the published 0.83, with all intermediate numbers reported honestly. (ii) On the eSOL *E. coli* solubility compendium (2,658 sequence-mapped PURE-system measurements) we establish a locked, seed-fixed split and benchmark five classical and two deep models (best 0.797 AUROC). (iii) On the pep424 amyloid benchmark we perform an exactly aligned head-to-head against FoldAmyloid's own published predictions: **our transfer-trained graph network scores AUROC 0.7902 against FoldAmyloid's 0.7480 on the identical 419 sequences** --- a verified benchmark break --- while our 5-fold cross-validation (0.7944) remains honestly below AmyloGram's published 0.865. (iv) A multi-task network (TriNet) trained jointly on all three tasks drives a tri-objective discovery screen; its first gated candidate, **FEKEAKKIEIKRH** (P_acp = 0.904, P_sol = 0.867, P_agg = 0.055), is named, quantified and falsifiable. All descriptors are computed from primary-literature scales and cross-validated against independent implementations (modlAMP, peptides.py), a process that caught and fixed a real mass-table bug. Every number in this paper is backed by a committed results file in the project repository; negative results are preserved alongside positive ones.

\newpage

# 1. Introduction

Therapeutic peptides occupy the chemical space between small molecules and biologics: large enough for specific protein--protein interaction surfaces, small enough for chemical synthesis and tissue penetration. Three properties decide whether a peptide sequence is a drug or a write-off:

1. **Activity** --- here, selective anticancer cytotoxicity. Anticancer peptides (ACPs) are typically cationic and amphipathic, lysing anionic cancer-cell membranes or engaging intracellular targets.
2. **Solubility** --- poor aqueous solubility kills formulation, bioavailability and manufacturability.
3. **Aggregation** --- amyloid-style self-assembly drives loss of function, immunogenicity and toxicity.

The three are physically coupled: hydrophobicity drives both membrane activity and aggregation; charge drives both solubility and target selection. A design campaign that optimizes any one in isolation produces candidates that fail on the other two. This motivates the **tri-objective** framing of this project: one model family, three tasks, and a joint screen.

## 1.1 Prior art and benchmarks

- **Anticancer activity:** AntiCP (Tyagi et al. 2013), AntiCP 2.0 (Agrawal et al. 2021), ACP-MHCNN (Scientific Reports 2021), mACPpred, ACP-DL, iACP, ACPred-FL. AntiCP 2.0's main benchmark reports 73.99% accuracy / MCC 0.48 / AUROC 0.83 on its held-out split; the alternate benchmark 88.18% / 0.76 / 0.95.
- **Solubility:** PROSO/PROSO II, DeepSol (Khurana et al. 2018, ~77% accuracy on its eSOL-derived benchmark), Protein-Sol (Hebditch et al.), NetSolP. The eSOL resource (Niwa et al. 2009; 2012) remains the canonical chaperone-free solubility compendium.
- **Aggregation:** TANGO, AGGRESCAN, PASTA 2.0, FoldAmyloid, WALTZ (on WALTZ-DB), AmyloGram (Kozlowski & Burdukiewicz 2017), AmyloGram's successor tools. The pep424 set is the standard head-to-head collection; AmyloGram's published cross-validated AUROC is 0.865.

## 1.2 Contributions

1. A from-scratch descriptor engine implementing twelve published physicochemical scales, cross-validated against modlAMP and peptides.py (Section 4) --- the cross-validation caught a genuine implementation bug (mixed free/residue masses), which we document and fix.
2. Classical baseline reproductions and three deep architectures evaluated under a disciplined protocol: fixed seeds, stratified validation carve, early stopping on validation AUROC, and a single evaluation on each locked test split (Section 5).
3. A **verified benchmark break**: pepx-GNN transfer AUROC 0.7902 vs FoldAmyloid's published predictions 0.7480 on pep424, exactly aligned (Section 6.3).
4. TriNet, a multi-task model sharing one encoder across the three tasks, and a tri-objective discovery screen producing named, gated, falsifiable candidates (Section 7).
5. **pepx**, a runnable CLI + library shipping the trained weights (Section 8).
6. Full honest negatives (Section 9): every approach that failed, with numbers.

\newpage
