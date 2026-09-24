# MEGA27-09c: Tri-Objective Peptide Intelligence (solubility / aggregation / anticancer)

Part of the MEGA-PROGRAM-27 computational-biology program (item 9, part 3).

## Goals
1. Benchmark-breaking classifiers for peptide anticancer activity (vs AntiCP 2.0 / ACP-MHCNN), solubility (vs DeepSol / Protein-Sol), and aggregation propensity (vs PASTA 2.0 / AGGRESCAN-class baselines), each on the published benchmark's own splits.
2. Discovery: a multi-task encoder + generative screen producing NAMED novel peptide candidates predicted simultaneously anticancer-active, soluble, and low-aggregation (tri-objective), quantified and falsifiable.
3. A reusable tool: `pepx` CLI + library for tri-objective peptide scoring.

## Rules of the house
- Real open datasets only; every dataset logged in data/MANIFEST.md with source URL + checksum.
- Hermetic pytest suite; no network at test time.
- CNN + GNN core architectures (PyTorch, CPU).
- All claims verified against locked splits; honest negatives preserved in results/.
