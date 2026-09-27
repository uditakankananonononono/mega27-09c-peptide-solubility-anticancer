# PREREGISTRATION C7 - Confirmatory constrained screen (locked 2026-09-27,
# BEFORE any confirmatory scoring)

Verdict item C7: "Pre-register the constrained screen BEFORE scoring
(thresholds, family, decision rules in this addendum series)."

The existing results/discovery_constrained.json (commit 1d13f33: 2,100
proposals, 35 candidates, 10 passing all gates) is hereby reclassified as
the EXPLORATORY PILOT batch - it was scored before any preregistration. It
stands in the record as pilot data only. The CONFIRMATORY batch below is
fresh-seeded and untouched at lock time. Bars below are pilot-informed;
that fact is disclosed, not hidden.

## Proposal family (locked; identical machinery to the pilot, disclosed)
- Model: shipped v0.1 TriNet, frozen weights results/trinet.pt +
  results/trinet_norm.npz. No retraining between pilot and confirmatory.
- Chain: single-point-mutation Metropolis, constrained kernel (descriptor
  gates reject BEFORE scoring): net_charge in [2, 9], gravy < 0, boman in
  [1.5, 3.2]. Seeds: the 60 CancerPPD sequences of length 13-16 that pass
  the descriptor gates (same seed pool as pilot); 60 steps per seed;
  temperature T = 0.05*(1 - step/90) + 0.005; acceptance on
  tri_score = p_acp * p_sol * (1 - p_agg).
- Confirmatory rng seed: 101 (pilot used 11). One batch, no re-seeding
  after seeing outcomes.

## Candidate gates (locked; carried over from the pilot unchanged)
- p_acp >= 0.9, p_sol >= 0.7, p_agg <= 0.2, plus the descriptor gates above.
- Novelty: sequence absent from CancerPPD + APD3 + AntiCP2 main + alternate
  (the pilot's known-pool rule; the DBAASP screen pool stays a separate
  check, and C6 BLAST novelty is reported separately per candidate).

## Decision rules (locked)
- SUCCESS: >= 10 candidates pass all gates in the confirmatory batch
  (matching the pilot's 10 - a replication bar, not an improvement bar).
- PARTIAL: 4-9 passing.
- FAILURE: <= 3 passing (reported as-is; negatives never terminal).

## Falsifiers / honesty controls (locked)
- Scramble control: every passing candidate gets one length-preserving
  shuffle (seed 103). The scramble pass rate should collapse toward the
  random panel rate; if >= 50% of scrambles still pass all gates, the gate
  set is flagged as composition-driven (same honesty pattern as 09a Q9).
- Random panel: 100 uniform-random 13-16-mers (seed 107) scored through the
  identical pipeline; expected pass rate ~0. If the random panel passes at
  or above the candidate rate, the screen result is void until diagnosed.

## Reporting
results/discovery_constrained_prereg.json + per-candidate table; pilot vs
confirmatory reported side by side; the pilot is never cited as
preregistered evidence.
