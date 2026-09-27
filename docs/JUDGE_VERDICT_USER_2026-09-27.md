# USER-PROVIDED JUDGE VERDICT - 2026-09-27 (user mega-verdict section 7)

Source (archive of record): the user's own WhatsApp message 2026-09-27 12:02:56
IST, wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMEFGRDY4MzY4OTkxNzFEQURGRAA=,
verified author=user (16,557 body bytes, 10-project mega-verdict). This repo's
section and the Cross-Cutting Computational Themes were extracted DIRECTLY from
that authenticated message. Note: a courier-compiled file received the same
minute did NOT byte-match her message; it was rejected as archive source.
Header directive (verbatim first words of her message): "IGNORE ABOUT ISEF
DELIVERABLES, IMPROVE PAGE COUNT" - page-length weaknesses and storyboard
items are superseded; papers grow with substantive content. Under her rule
("EACH PROJECTS NEED ONE FROM ME TO PASS", WhatsApp 10:01:47 same day) this is
this project's ONE counted judge round. Factual claims are verified
independently before adoption.

---

7. paper_with_updates.pdf: Tri-Objective Peptide (09c)
Weaknesses (20) — computational only:

AntiCP 2.0 main: 0.8029 vs published 0.83 — not beaten.

AntiCP 2.0 alternate: 0.931 vs 0.95 — not beaten.

AmyloGram CV 0.865 stands; ours 0.7947.

Only verified win is FoldAmyloid on pep424 (+9.1 AUROC).

Discovery candidates unvalidated.

25-epoch screen drifted anionic — model failure.

Random-weight CNN falsifier not run.

Mass-table bugs (P/Q/V/W/Y scramble) survived test suite.

Shipped weights almost overwritten silently — integrity incident.

Single-seed headline numbers.

AntiCP 2.0 leakage: 16.6% test peptides near-homolog to train.

Composition descriptors dominate AntiCP benchmark.

CNN loses to 3-mer LR on CPP (0.8956 vs 0.9265).

GNN underperforms CNN on AMP tasks.

Tool count 24, not 40 (gate open).

Only 287 accession records, not 120 cohorts.

Limited to 10-30-mer screens.

Conformational blindness — no membrane state.

Paper is 47-52 pages — too long.

Constrained-screen family not pre-registered.

Additions (computational):

Add a hemolysis classifier (HemoPI-class) as a fifth task.

Add ESM-2/ProtT5 embedding baseline.

Add a disulfide-aware gate.

Add per-residue solubility profile (CamSol-style).

Raise pytest coverage from 29.2% to 60%+.

Add NCBI BLAST novelty check.

Pre-register the constrained screen.

Simplify: lead with pep424 break + constrained-screen method.

Reduce to 12 slides.

---

Cross-Cutting Computational Themes
Recurring weaknesses:

Dataset/tool count inflation (nested records counted as independent).

Long papers (47-58 pages) — not ISEF-ready.

Negative-heavy narratives that obscure positive contributions.

Single-seed headline numbers.

Ad hoc gates/thresholds rather than theoretically derived.

Homology leakage in random-split benchmarks.

No leave-family-out CV in most projects.

CIs often overlapping — point-estimate wins only.

Universal computational additions:

One primary question per paper.

One locked primary endpoint.

Cluster-level bootstrap CIs everywhere.

Leave-family-out CV as the primary protocol.

12-slide storyboard as the ISEF deliverable.

One-page summary card.

Decision tree for tool use.

Pre-registered replication within the paper itself.

"What this is NOT" section in every abstract.

Reduce tool/dataset counting to study-level units.
