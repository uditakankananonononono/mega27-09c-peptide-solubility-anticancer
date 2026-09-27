# PREREGISTRATION C3 - Disulfide-aware gate (locked 2026-09-27, BEFORE any
# C3 scoring)

Verdict item C3: "Disulfide-aware gate in the screen."

## Gate rule (locked)
- The constrained screen kernel gains one rule: a proposal's Cys count must
  be EVEN (0, 2, 4, ...). Rationale (disclosed as heuristic, not a
  structural claim): in disulfide-stabilized peptide scaffolds an odd Cys
  count leaves a reactive free thiol - a covalent/aggregation liability. The
  gate applies BOTH in-kernel (proposals rejected before scoring) and in the
  final candidate gate set.
- No other C7 machinery changes: same frozen TriNet v0.1, same descriptor
  gates (net_charge [2,9], gravy < 0, boman [1.5,3.2]), same candidate
  thresholds (p_acp >= 0.9, p_sol >= 0.7, p_agg <= 0.2).

## Batch (locked)
- Fresh rng seed 109 (C7 pilot used 11, confirmatory 101). Same 60-seed /
  60-step protocol, same seed pool.

## Decision rules (locked)
- SUCCESS: >= 5 candidates pass all gates including even-Cys (the stricter
  kernel still produces candidates).
- Otherwise reported as-is (PARTIAL 1-4, NONE 0). This is a screen
  refinement analysis, not a benchmark beat; the C7 confirmatory result
  stands untouched.
- Overlap with the C7 confirmatory top-10 reported (expected low - fresh
  seed, stricter kernel).

## Falsifiers (locked)
- Scramble control seed 111 and random panel seed 113, same construction and
  interpretation as C7; the even-Cys rule applies to scrambles/panel too.

## Reporting
results/c3_disulfide_screen.json, as-is, side by side with the C7 numbers.
