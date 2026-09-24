# 3. Data

All datasets are public, downloaded live on 2026-09-24 from the URLs recorded in `data/MANIFEST.md` (repository root). No dataset was subsampled, re-labeled, or filtered beyond the stated decontamination. Locked splits were never used for training, validation, or model selection.

## 3.1 Anticancer peptides

**AntiCP 2.0 (Agrawal et al. 2021).** Downloaded from the authors' server (`webs.iiitd.edu.in/raghava/anticp2/download.php`). Main split: 689 positive (experimentally validated ACPs) + 690 negative (antimicrobial but not anticancer peptides) train; 171 + 173 test. Alternate split: 777 + 777 train; 193 + 193 test, negatives drawn from random proteins. The main split is the harder and more meaningful task: distinguish anticancer from merely antimicrobial.

**CancerPPD (Tyagi et al. 2015).** 2,647 canonical L-amino-acid ACPs (`l_natural.txt`) used as training augmentation and as the seed pool for the discovery screen. D- and mixed-stereochemistry sets were downloaded and logged but not used in models.

**APD3 (2024 release).** 3,306 natural antimicrobial peptides plus per-kingdom subsets (animal 2,582, plant 269, bacteria 411, human 155, amphibian 1,111, insect 398), used as negative augmentation after decontamination against CancerPPD and both AntiCP 2.0 splits.

## 3.2 Solubility

**eSOL (Niwa et al. 2009/2012, LSDB archive).** 4,132 *E. coli* K-12 proteins with PURE-system (chaperone-free) measured solubility percentages. Sequences were mapped from UniProt proteome UP000000625 by gene name, yielding 2,658 sequence-labeled records; binary labels use the published 30% cutoff (62.3% soluble). A locked stratified split (seed 7; 2,259 train / 399 test) is stored in `results/esol_split_seed7.json`.

## 3.3 Aggregation

**AmyloGram / WALTZ-DB-derived sets.** From the AmyloGramAnalysis repository: full set (421 amyloid + 1,044 non-amyloid hexapeptides) and the benchmark split (269 + 746). Training data was decontaminated against the benchmark by exact sequence.

**pep424.** 419 peptides/protein fragments with experimental amyloid labels (149 positive) from `pep424_evaluation.txt`. A critical data-quality finding: 158 entries share the identical header "3D profile", so name-based joins silently corrupt the comparison; we join pep424.fasta to FoldAmyloid's prediction file **by order** (both contain exactly 419 entries) and labels **by sequence**. This is documented in Section 6.3 and was the difference between an invalid and a valid head-to-head.

**FoldAmyloid and PASTA 2.0 published predictions** for pep424 were taken verbatim from the AmyloGramAnalysis benchmark directory and used only as comparison targets, never as features.

## 3.4 Decontamination protocol

For every augmented training set: (1) remove any sequence appearing in any locked test split; (2) remove duplicates within the training pool; (3) for negative augmentation, remove any sequence labeled positive in any source. Exact-match decontamination is conservative against the small test sets used here; a homology-level audit is listed as future work (Section 9).

\newpage

# 4. Descriptor engine and cross-validation

The descriptor engine (`pepx.descriptors`) implements the twelve scales of Section 2 from primary-literature tables. Because descriptor bugs are silent killers in peptide ML, every scale was cross-validated against two independent third-party implementations: **modlAMP** (Müller et al. 2017) and **peptides.py**.

## 4.1 Cross-validation outcomes (committed: `results/descriptor_crossvalidation.json`)

| quantity | pepx | independent | verdict |
|---|---|---|---|
| Instability index | -14.97 / 84.68 / 44.73 / 294.0 | modlAMP: identical to 2 dp | exact match |
| Molecular weight | (after fix) 1524.0 | modlAMP 1524.0 | match after bugfix |
| Isoelectric point | 10.70 | modlAMP 13.23 | pKa-set dependence, documented |
| GRAVY | 0.586 | peptides.py -0.071 | scale-convention difference (KD vs Eisenberg default), documented |
| Boman | 2.561 | peptides.py 0.456 | normalization-convention difference, documented |

**Bug found and fixed.** The first cross-validation run disagreed with modlAMP on molecular weight by ~90 Da on a 14-mer. Inspection showed our amino-acid mass table mixed free-amino-acid and residue masses for 11 of 20 residues. The corrected table (free average masses, verified residue-by-residue) reproduces modlAMP exactly; a regression test (`test_molecular_weight_matches_modlamp`) locks the fix. This episode is exactly why the cross-validation protocol exists, and it is reported here rather than silently patched.

**pI is pKa-set-dependent.** modlAMP/EMBOSS pKa values place KLAKLAKKLAKLAK at pI 13.23; the Bjellqvist set places it at 10.70. Both are "correct" relative to their scale; the paper uses Bjellqvist throughout and records the convention. Similarly GRAVY (Kyte--Doolittle) and peptides.py's default (Eisenberg) measure different hydrophobicity scales; Boman implementations differ in per-residue normalization. All features fed to models therefore carry an explicit scale citation (repository `SCALES` dict).

\newpage
