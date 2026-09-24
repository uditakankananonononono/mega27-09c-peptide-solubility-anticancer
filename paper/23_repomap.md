# Appendix I. Map of the repository

```
mega27-09c-peptide-solubility-anticancer/
  src/pepx/            the library: alphabet, fasta, descriptors (12 published
                       scales + dual-validated mass table), encoders (AAC/DPC/
                       k-mer/descriptor), datasets (5 loaders, strict
                       validators), models (CNN/GNN/CNNv2/TriNet), trainer
                       (locked-seed disciplined loop), cli, api (FastAPI)
  tests/               29 hermetic tests: contracts, identities, both mass
                       locks, property-based suites, CLI end-to-end
  runs/                every experiment script, in execution order; logs kept
                       alongside (run_*.log), including the failed ones
  results/             every number in this paper as committed JSON, plus the
                       shipped artifacts: trinet.pt (v0.1 weights),
                       trinet_norm.npz (train statistics), trinet.onnx
                       (portable export), trinet_25ep.pt (drifted negative)
  data/                raw/ (benchmark files, untouched) + MANIFEST.md (live
                       source URLs, fetch dates, accession lists) +
                       uniprot_accessions/ (128 entries) + external/ (DBAASP,
                       ChEMBL, PDB, attempt logs)
  paper/               this document's sources (00-23, pandoc -> paper.pdf)
                       and figures fig1-fig8
```

Reading order for a new contributor: `data/MANIFEST.md` (what the data is), `src/pepx/descriptors.py` (the scales and their citations), `runs/run_pep424_v3.py` (the break), `runs/run_discovery_constrained.py` (the screen), `src/pepx/cli.py` (the tool), then this paper front to back. The git history is append-only in spirit: no results file was ever deleted, and the two history rewrites on record (commit-message page-count corrections) are disclosed here.
