# Roadmap

Ordered by expected value per unit effort, each item with its acceptance test.

1. **Hemolysis model (HemoPI-class).** The candidate family is cationic; hemotoxicity is the binding safety gap. Acceptance: locked-split AUROC reported against HemoPI's published numbers, same discipline as Section 5. Data: Hemolytik (download endpoint currently non-machine-readable - will parse the per-entry pages or fall back to the HemoPI training files).
2. **ESM-2/ProtT5 embedding baseline.** Expected +1-3 AUROC on all three tasks (published on comparable tasks). Acceptance: same locked splits, same seed, honest win or loss recorded. Requires more compute than this container; queued for any GPU day.
3. **Disulfide-aware gate.** Two of the ten constrained-screen leaders carry cysteine; oxidative dimerization is unpriced. Acceptance: gate added, candidate family re-screened, deltas reported.
4. **Per-residue solubility profile (CamSol-style).** Turns the scalar P_sol into a profile for lead optimization. Acceptance: profile correlation vs CamSol on the accession panel.
5. **Coverage to 60%.** Current 29.2% (honest). Target: trainer and model paths under property tests. Acceptance: pytest-cov report committed per release.
6. **NCBI BLAST novelty check.** The 9,683-sequence screen is a pool check, not a proteome check. Acceptance: nr/swissprot BLAST of the candidate family, identity table committed.

# Glossary and notation

**AAC** amino-acid composition (20 features). **DPC** dipeptide composition (400). **AUROC** area under the receiver operating curve; equals the probability a random positive outranks a random negative (Section 19.1). **MCC** Matthews correlation coefficient; balanced binary quality at a fixed threshold. **GRAVY** grand average of hydropathy (Kyte-Doolittle means over the sequence). **Boman index** estimated protein-binding potential (kcal/mol scale, sign conventions per Section 2). **pI** isoelectric point (Bjellqvist pKa set here, bisection root of net charge, proof in Section 2). **DIWV** Guruprasad instability dipeptide weights. **Locked split** a train/test partition fixed by seed and stored (`results/esol_split_seed7.json`), never re-randomized. **Aligned comparison** a head-to-head where both models' scores cover the identical sequences in the identical order (pep424 protocol). **Gate** a descriptor-range constraint applied at candidate acceptance (charge, GRAVY, Boman, length) or, in the constrained kernel, at proposal time. **Accession-level dataset** one identifier-backed record individually fetched and used (program counting rule). **Negative result** a completed, committed experiment whose outcome falsifies or bounds an approach; in this repo negatives are kept, named and cited, never deleted.
