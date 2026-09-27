# PREREGISTRATION ADDENDUM - 2026-09-27 (user mega-verdict, section 7)

Source: user WhatsApp 2026-09-27 12:02:56 IST (wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMEFGRDY4MzY4OTkxNzFEQURGRAA=, verified author=user). Verbatim archive: docs/JUDGE_VERDICT_USER_2026-09-27.md.
Queue LOCKS at commit.

## Locked queue (compute, cheap-first)
- C1 Hemolysis classifier (HemoPI-class dataset) as a fifth task.
- C2 ESM-2/ProtT5 embedding baseline on existing splits.
- C3 Disulfide-aware gate in the screen.
- C4 CamSol-style per-residue solubility profile (in-house implementation
  against published method description; disclosed as re-implementation).
- C5 Pytest coverage 29.2% -> 60%+ (locked target).
- C6 NCBI BLAST novelty check for discovery candidates (blastp short-seq vs
  nr; if egress blocked, use local Swiss-Prot blast - fallback disclosed).
- C7 Pre-register the constrained screen BEFORE scoring (thresholds, family,
  decision rules in this addendum series).
- C8 Random-weight CNN falsifier (explicitly missing per verdict).
- C9 Paper leads with pep424 break + constrained-screen method.
- In-flight DeepSol arm (predates this verdict) continues and reports as-is.
- Page-length/storyboard weaknesses SUPERSEDED by her header directive.

## Factual verification notes
- AntiCP2 0.8029/0.931 vs published, AmyloGram 0.865 vs ours 0.7947,
  FoldAmyloid +9.1 AUROC on pep424: consistent with committed results.
- 25-epoch anionic drift, mass-table scramble, near-overwrite incident:
  confirmed in repo history; already documented as failures/incidents.
- "Tool count 24 not 40", "287 accession records not 120 cohorts": recount at
  study-level units locked (cross-cutting); paper corrected.
