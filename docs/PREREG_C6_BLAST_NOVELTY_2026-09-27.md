# PREREGISTRATION C6 - BLAST novelty check (locked 2026-09-27, BEFORE any C6
# query)

Verdict item C6: "NCBI BLAST novelty check for discovery candidates (blastp
short-seq vs nr; if egress blocked, use local Swiss-Prot blast - fallback
disclosed)." Live egress test 2026-09-27: NCBI BLAST URL API reachable (200),
so the PRIMARY path is NCBI blastp vs nr; the local Swiss-Prot + mmseqs2 arm
runs as a disclosed cross-check, not a fallback.

## Candidate set (locked)
- The C7 confirmatory batch's 10 passing candidates (preregistered, 64c88d7).
- The C3 disulfide-gated batch's 12 passing candidates (preregistered, bbd590e).
- The C7 exploratory PILOT batch's 10 passing candidates (disclosed as pilot,
  reported separately).
No additions or removals after results.

## Query protocol (locked)
- NCBI URL API: PROGRAM=blastp, DATABASE=nr, one RID per candidate, >=10 s
  between submissions (NCBI usage guidelines), polled with >=30 s spacing.
  If NCBI rejects or truncates short-sequence queries, that is reported
  as-is (the Swiss-Prot arm still stands).
- Cross-check: UniProtKB/Swiss-Prot (current release fasta, sha256 recorded),
  mmseqs2 easy-search, sensitivity -s 7.5, alignment mode covering full
  query length, e-value <= 10.

## Novelty rule (locked, applied identically to both arms)
- A candidate is NOVEL if it has no hit with >= 90% sequence identity over
  >= 80% of its length (any species/entry). Near-misses (best hit identity
  70-90%) are reported as NEAR with the top hit's accession, identity, and
  coverage. Hits >= 90%/80% = NOT NOVEL, reported with accessions.
- Low-complexity cationic peptides are expected to have many partial hits;
  the identity-over-full-length rule, not hit count, decides novelty.

## Decision (locked)
- Per candidate: NOVEL / NEAR / NOT NOVEL with evidence. Screen-level: the
  fraction NOVEL is reported as-is; there is no beat bar (this is a novelty
  audit, not a benchmark).
