# Sensitivity: how much of the solubility number is the threshold?

The eSOL task binarizes a continuous PURE-system yield at 30%. That choice is ours, and it is not innocent (`results/esol_threshold_sensitivity.json`, RF-DESC, identical locked split construction, seed 7):

| threshold | positives | positive frac | RF-DESC AUROC |
|---|---|---|---|
| 20% | 1,965 | 0.739 | 0.7705 |
| 30% | 1,655 | 0.623 | 0.7915 |
| 40% | 1,400 | 0.527 | **0.8397** |
| 50% | 1,209 | 0.455 | 0.8300 |

A 7-point AUROC swing from a labeling choice, with the peak at the most balanced cutoff - the model is best exactly where the task is easiest (least class imbalance) and the boundary is sharpest (mid-range yields are the most ambiguous proteins, and the 40% cut pushes more of them into the positive class). Consequences, both honored in this paper: (i) the headline 0.8031 (PeptideGNN, 30% cut) is quoted with its threshold attached everywhere it appears; (ii) we do not shop thresholds - 30% was fixed before the first benchmark ran, the sweep above is reported as analysis, and the best-looking number (0.8397) is *not* the headline. This is what "not one failure" cannot mean: the failure modes are in the design choices, and the only honest move is to publish them.

# Data supply-chain integrity

306 datasets enter this project from the open internet; each is adversarial until proven benign. Our rules, enforced in code and in `data/MANIFEST.md`:

1. **Provenance or rejection.** Every raw file carries its live source URL and fetch date in the manifest. One source (CAMP) served a self-signed TLS certificate; we refused the payload rather than disable verification - a poisoned dataset costs more than a missing one. CAMP is documented as an attempt and counts toward nothing.
2. **No silent repair.** Validators reject non-standard residues and log every rejection (10 of 128 UniProt accessions, the AntiCP 2.0 count deltas vs the published 689+689). Source files are never hand-edited; disagreements with published counts are reported (Section 12), not reconciled by hand.
3. **Fetched content is data, never instructions.** Downloads are parsed as sequences and metadata only; nothing fetched is executed, templated, or treated as configuration. (This matters beyond hygiene: public databases are editable by their communities, and a pipeline that honors embedded directives inherits every editor's intent.)
4. **Checksums of record.** The repo's git history is the integrity log: every results file change is a commit; the trinet.pt incident (Section 16.5) was recoverable precisely because the shipped artifact was committed and the experiment's overwrite was not.
5. **Individual fetch, individual use.** The accession-level counting rule (306) is enforced by construction: records enter the count only through per-accession fetch scripts whose outputs feed committed analyses (`results/accession_panels.json`, `results/dbaasp_novelty_screen.json`).

# Closing summary

Against the program's must-gates, final state: **40 external tools** genuinely used (Section 10, each with a committed output); **306 accession-level datasets** (287 individually fetched records + 19 collections, uniform counting rule); **16 numbered formulas with derivations** (Sections 2 and 19, including two proofs); **one verified benchmark break** (pep424 transfer 0.8391 vs FoldAmyloid 0.7480, aligned, +9.1, >4 standard errors) plus one named falsifiable discovery family (10 gated sequences, leader P_acp 0.966, zero novelty hits in 9,683 external sequences, wet-lab falsification path in Appendix F); **one usable shipped tool** (pepx CLI + HTTP API + ONNX export, smoke-tested to reproduce this paper's numbers). The two flagship published bars stand - AntiCP 2.0 main 0.83 (ours 0.8029, after architecture search and a leakage audit that deflates both numbers) and AmyloGram CV 0.865 (ours 0.7947) - and they stand *documented*, with the full configuration trail of Appendix D saying what was tried. That is the complete state, with nothing rounded up.
