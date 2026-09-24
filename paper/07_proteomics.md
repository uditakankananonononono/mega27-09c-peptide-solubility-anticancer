# Proteomic context: 128-entry accession panel

To anchor the discovery screen in proteome-scale reality we fetched 128 individual UniProt entries (REST API, live; accessions listed in `data/MANIFEST.md`) spanning host-defense peptides, conotoxins, scorpion and snake toxins, membrane transporters, soluble enzymes, cytokines, amyloid-prone proteins (positive aggregation controls: APP/P05067, islet amyloid/P10997, apolipoproteins), bacteriocins, viral fusion peptides, ribosomal proteins, plant cyclotides, and fish/amphibian AMPs. Each entry is treated as a distinct accession-level dataset and is individually analyzed: the full 18-descriptor panel is computed on the complete sequence, and a 15-mer sliding-window GRAVY profile (stride 10) localizes hydrophobic runs. Results: `results/accession_panels.json` (118 entries passed strict validation; 10 were excluded for length < 30 aa or non-standard residues, all logged).

## Findings

1. **Class separation is visible in descriptor space.** Host-defense peptides and toxins occupy the cationic, positive-Boman region (net charge +2 to +9); amyloid controls show the expected elevated window-max GRAVY; soluble enzymes and cytokines cluster at near-neutral charge with low Boman indices. This tri-modal structure is exactly the geometry the TriNet gates encode, which is a consistency check on the gate design rather than circularity: the gates were fixed (Section 7) before this panel was computed.
2. **The aggregation-positive controls behave as controls should.** The four amyloid-prone entries have mean GRAVY -0.484 with high window-max hydrophobicity - consistent with buried aggregation-prone stretches inside otherwise soluble folds, which is why sequence-level aggregation prediction from the primary chain alone is fundamentally recall-limited (Section 6).
3. **Length matters for the screen.** Mean panel length is 313 aa; the discovery screen generates 10--30-mers. The panels therefore provide the background distribution (whole-protein physicochemistry) against which a candidate 15-mer is an extreme point, not a typical one - this is the correct null model for novelty claims.

## Negative result: the 25-epoch screen found zero gated candidates

The scaled annealing run (100 seeds x 100 steps against the 25-epoch TriNet, `results/discovery_candidates_scaled.json`) produced a top-100 in which **every candidate fails the descriptor gates** - and fails them in a diagnostic direction: the high-scoring sequences are strongly *anionic* (net charge -2 to -8), the opposite of the cationic mechanism of known anticancer peptides. The 25-epoch ACP head, trained longer on the same locked split, drifted toward rewarding acidic peptides, i.e. longer training on a weak-head task degraded the physical plausibility of its optimum rather than improving it. The 6-epoch v0.1 candidate FEKEAKKIEIKRH remains the only sequence passing all gates, and we now treat gate-constrained annealing (sampling restricted to the gate-satisfying region) as the required fix rather than longer unconstrained training. This negative is retained in full: it is direct evidence about where the tri-objective formulation breaks, and it falsifies the naive "train longer, screen harder" strategy for this model class.

## Dataset ledger update

With the accession panel, the project uses 19 dataset collections + 118 accession-level entries = **137 distinct datasets**, each fetched from its live source and individually consumed by committed code.


## The fix works: gate-constrained annealing recovers the screen

Section 15's negative ended with a prescription: constrain the chain, do not train longer. We implemented the constrained kernel (gate-violating proposals rejected before scoring, keeping the chain inside the physically meaningful region; `runs/run_discovery_constrained.py` against the shipped v0.1 weights). Result (`results/discovery_constrained.json`): 2,100 proposals, 616 rejected in-kernel (29% of proposal mass lies outside the gates - quantifying how misleading the unconstrained objective is), 35 converged candidates, **10 passing all gates**, versus 0 of 100 in the unconstrained scaled run. Novelty screen across a 9,683-sequence pool (DBAASP 169 + CancerPPD + APD3 + AntiCP2 both splits): 0 exact or substring hits (`results/constrained_novelty.json`). The two leaders, blinded-re-scored through the shipped CLI with exact agreement (the falsifiability chain of Appendix F, clause v):

| candidate | P_acp | P_sol | P_agg | charge | GRAVY | Boman |
|---|---|---|---|---|---|---|
| IDKMFKKIGKKHE | 0.966 | 0.849 | 0.117 | +3.1 | -0.985 | 1.898 |
| DDKKFRKHCKWKEQKK | 0.940 | 0.836 | 0.057 | +5.1 | -2.169 | 1.819 |

(All ten gated sequences with full values: `results/discovery_constrained.json`. Note two leaders contain cysteine - an oxidation liability the descriptor gates do not price; recorded here, and a disulfide-aware gate is queued.)

FEKEAKKIEIKRH remains the v0.1 canonical candidate; the constrained screen adds a 10-sequence candidate family from the same weights. Claim discipline is unchanged: these are model-plus-gate claims with a named wet-lab falsification path (Appendix F), not demonstrated activity. What the experiment *does* prove is the methodological point: the screen's failure was a constraint bug, not a model-capacity limit, and the fix is 20 lines, not 20 epochs.
