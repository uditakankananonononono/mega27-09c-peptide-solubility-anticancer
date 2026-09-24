# Using the tool: a practical guide

## Install and first score

```
git clone <repo>; cd mega27-09c-peptide-solubility-anticancer
pip install torch numpy scikit-learn scipy pandas biopython   # core deps
PYTHONPATH=src python3 -m pepx.cli score FEKEAKKIEIKRH
# sequence  P_anticancer  P_soluble  P_aggregation
# FEKEAKKIEIKRH  0.9033  0.8667  0.0550
```

Batch mode: `score-file input.fasta --out scores.tsv` (~2,000 sequences/min). Descriptors only: `describe SEQ` emits the full 18-value panel as JSON. HTTP: `uvicorn pepx.api:app` then POST `/score` with `{"sequences": [...]}` (batch limit 512; validation rejects non-standard residues with a named-letter diagnostic). ONNX: `results/trinet.onnx` runs in onnxruntime without PyTorch (parity 4.8e-7 logits; feed train-standardized descriptors, mu/sd in `results/trinet_norm.npz`).

## Reading the numbers

The three probabilities come from three different label universes (Section 12) and are *not* calibrated against each other: P_acp 0.9 does not mean "90% chance of killing cancer cells", it means "deeper into the ACP region of training space than 90%-equivalent of the calibration-free score range". The descriptor gates are the semantic layer: a sequence is only a *candidate* when it also passes charge +2..+9, GRAVY < 0, Boman 1.5-3.2, length 10-30. The shipped smoke test (API/CLI agreement on FEKEAKKIEIKRH and KLAKLAKLAKLAK) is the regression tripwire - if your install disagrees with those four numbers, the weights/normalization pair is mismatched, as happened during the Section 16.5 incident.

## Failure modes a user will actually meet

1. **Sequences with U/X/B/Z/J** are rejected, not silently corrected; map them first.
2. **Long proteins** (>60 aa) are truncated by the shipped model's input width; score windows instead (the accession-panel script shows the 15-mer stride-10 pattern).
3. **Cysteine- and methionine-rich candidates** carry unpriced liabilities (Section 15.1); treat their scores as optimistic.
4. **Scores far outside the training envelope** (very long, highly repetitive, all-aromatic sequences) are extrapolations; the leakage audit (Section 18) shows how much even *in-distribution* scores rely on near-homology.

# Statistical discipline

Every inferential number in this paper follows four rules, stated here so the reader can audit us:

**Rule 1 - the test set is a singleton event.** Locked test data is evaluated once per run, after all selection (early stopping, HPO, ensembling) is complete. Anything else - including choosing the best of 12 HPO trials by test score - is selection on test, and Section 5.4 quantifies the bias it would have created (0.731-0.788 spread, i.e. up to +5.7 points of fake "progress").

**Rule 2 - intervals travel with point estimates.** AUROC intervals use DeLong variance (Section 19.1 derivation); proportion intervals use Wilson scores (statsmodels, `results/tool_analyses.json`). Where we report a difference (the +9.1 pep424 break), the interval arithmetic is done on the difference, not eyeballed from the two point estimates.

**Rule 3 - multiplicity is acknowledged, not Bonferroni'd into hiding.** We ran ~19 model configurations (Appendix D). With that many looks, a nominal p < 0.05 "win" would be expected by chance; our break claim therefore rests on effect size (+9.1 AUROC, > 4 standard errors) rather than a thresholded p-value, and our near-misses (AntiCP2 0.80 vs 0.83) are reported as losses, not as "trending".

**Rule 4 - single-seed means single-seed.** Headline numbers are one seed (7). Where variance was measured it is printed next to the claim (10-config sweep, 12-trial HPO). Where it was not, the claim's language is bounded ("reaches", never "achieves robustly") and Limitations says so.

# Descriptor scale provenance

The 12 computed descriptor families and their primary-literature sources, as declared in the repository `SCALES` registry (verbatim):

| scale | source |
|---|---|
| KD_HYDROPATHY | Kyte & Doolittle, J Mol Biol 157:105-132 (1982) |
| EISENBERG | Eisenberg et al., Faraday Symp Chem Soc 17:109 (1982) |
| HOPP_WOODS | Hopp & Woods, PNAS 78:3824 (1981) |
| CF_HELIX / CF_SHEET | Chou & Fasman, Adv Enzymol 47:45-148 (1978) |
| AGGREGATION_ZZ | Zhou & Zhou, FEBS Lett 557:21-26 (2004) |
| DIWV (instability) | Guruprasad, Reddy & Pandit, Protein Eng 4:155-161 (1990) |
| BOMAN | Boman, Antimicrob Agents Chemother 47:1904 (2003) |
| ALIPHATIC | Ikai, J Biochem 88:1895-1898 (1980) |

plus derived quantities with first-principles definitions in Section 2: net charge at pH 7 (Henderson-Hasselbalch sums), pI (bisection root, uniqueness proof), molecular weight (free-average masses minus water, dual-validated), aromaticity/fractions (counting identities), isoelectric and aliphatic indices. Every scale in a model input names its provenance; the two-bug saga of Section 4 is what happens when one silently does not.
