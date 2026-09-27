# PREREGISTRATION C4 - CamSol-style per-residue profile (locked 2026-09-27,
# BEFORE any C4 fitting or scoring)

Verdict item C4: "CamSol-style per-residue solubility profile (in-house
implementation against published method description; disclosed as
re-implementation)."

## Method (functional form from the published description)
- Per-residue linear combination s_i = a_H*p_H(i) + a_C*p_C(i) + a_a*p_a(i)
  + a_b*p_b(i), then 7-residue window mean smoothing; overall intrinsic
  score = mean of the smoothed profile; sign convention: negative =
  insoluble (Sormanni et al., JMB 2015, doi:10.1016/j.jmb.2014.09.026,
  main text retrieved 2026-09-27, www-vendruscolo.ch.cam.ac.uk PDF; method
  re-stated in CSH Perspect Biol 2019 a033845).
- Scope cut (locked, disclosed): the hydrophobic-pattern and gatekeeper
  correction terms are described in the paper but NOT implemented in this
  v1 re-implementation.
- The paper's pattern/gatekeeper equations render as images in accessible
  copies; the numeric weights (Table S3) live in the paywalled JMB
  supplementary, which this sandbox cannot reach (fetch attempts logged).
  DISCLOSED: the 4 linear weights are therefore FIT IN-HOUSE on the eSOL
  seed-7 TRAINING split only (ridge, closed form). This artifact is a
  CamSol-STYLE re-implementation, never labeled CamSol.

## Scales (all published constants, cited in src/pepx/descriptors.py SCALES)
- p_H: Kyte-Doolittle hydropathy (1982), z-normalized across the 20 AA.
- p_C: charge at neutral pH (D,E = -1; K,R = +1; H = +0.1; else 0).
- p_a, p_b: Chou-Fasman helix/sheet propensities (1978), z-normalized.

## Validation (locked)
- Split: the committed seed-7 eSOL split (results/esol_split_seed7.json).
- Primary: Spearman correlation of the overall score vs measured
  solubility % on the TEST split.
- Secondary: AUC at the published 30% solubility cutoff vs the repo's
  committed RF-descriptor baseline AUC 0.7915 (same split,
  esol_threshold_sensitivity.json). BEAT: AUC > 0.7915 AND 10k bootstrap
  lower95 > 0.7915; MATCH: CI spans; MISS: below. Reported as-is.
- Falsifier: weights fit on seed-13 shuffled training labels must give
  test |Spearman| < 0.1, else the harness is void until fixed.

## Reporting
results/c4_camsol_profile.json (weights, metrics, CIs, falsifier) +
per-residue profiles for two disclosed example sequences.
