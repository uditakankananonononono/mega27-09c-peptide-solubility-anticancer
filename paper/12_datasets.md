# Dataset deep dives

Every collection used, with its provenance, its exact composition as loaded by `pepx.datasets` (post-validation counts; `results/dataset_stats.json`), and the role it plays. Numbers here are computed from the committed raw files, not quoted from source papers.

## AntiCP 2.0 (Agrawal et al. 2021, IIIT-Delhi)

Role: the primary anticancer benchmark. Main split as loaded: 860 positive / 863 negative (median lengths 18 and 26 aa, ranges 3-50 and 2-50). Alternate split: 969/972 (medians 19/28). The small deviations from the published 689+689/172+172 come from our strict validator (non-standard residues rejected and logged) and the published files containing a handful of duplicates and non-standard letters. Main-split negatives are experimentally characterized non-ACP antimicrobial peptides; alternate negatives are random UniProt-derived. That difference in negative provenance is why the alternate split scores higher for every model we ran - "does not look like a random protein fragment" is an easier decision boundary than "is an AMP but does not kill cancer cells" - and why we treat the main split as the honest bar.

## CancerPPD (Tyagi et al. 2015)

Role: positive-pool enrichment for discovery screening and logo analysis. 2,647 experimentally tested anticancer peptides after length filtering (median 16 aa, range 5-63). No labels beyond "tested ACP"; used only where a label-free positive pool is correct (screen background, k-mer statistics, cross-validation sources).

## eSOL / PURE-system solubility (Niwa et al. 2009; LSDB archive)

Role: the solubility task. 2,658 sequence-mapped proteins after our b-number mapping fix (Section 15): 1,655 soluble / 1,003 low-solubility at the 30% threshold, median lengths 243/317 aa. These are full proteins, not peptides - the model must learn solubility signals that transfer down to 10-30-mers, and the accession panel (Section 15) is the bridge evidence that the descriptor geometry is length-consistent.

## AmyloGram / WALTZ-DB (Burdukiewicz et al. 2020)

Role: the aggregation task and the unbeaten bar. 420 positive / 1,044 negative hexapeptide-centric entries (median length 6 both classes). Extreme class imbalance (1:2.5) and extreme shortness make this a different problem shape from the others; our GNN's 5-fold CV of 0.7947 vs the published 0.865 is reported unbeaten.

## pep424 + FoldAmyloid predictions

Role: the head-to-head benchmark. 424 hexapeptides with WALTZ-DB labels (419 cleanly aligned, 149 positive) and FoldAmyloid's own published prediction scores, order-aligned (Section 6.3; the 158-way "3D profile" header collision documented in Appendix G). Our transfer AUROC 0.8391 vs 0.7480 is the project's verified break.

## APD3 (Wang et al. 2016)

Role: augmentation negative/general-AMP pool. 3,203 natural AMPs (median 28 aa). Used for the augmentation experiment (Section 5, neutral result) and as a screen background.

## UniProt accession panel (128 fetched, 118 used)

Role: proteomic context. Individually fetched entries spanning 12 functional classes (Section 15), median length 228 aa; each individually analyzed (`results/accession_panels.json`). Ten exclusions (length < 30 or non-standard residues) are logged in the fetch script output.

## DBAASP (Pirtskhalava et al. 2021)

Role: external novelty screen. 169 individually fetched monomer records (median 15 aa); all 21 screened candidates show 0 exact or substring hits (`results/dbaasp_novelty_screen.json`).

## External evidence collections

RCSB PDB entries 2K6O and 2MAG (structure metrics, Section 17); ChEMBL (8 molecules, 60 activity records, potency table); AmyloGramAnalysis archive (benchmark files + competitor predictions); Biopython ProtParamData (DIWV instability table). Accession-level lower bound: 287 individually fetched and used records (118 UniProt + 169 DBAASP). The asserted 19 separate collections need per-collection verification before addition; see `DATASET_LEDGER.md`.
