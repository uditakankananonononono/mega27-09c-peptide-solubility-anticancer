# Paper-build evidence log, 2026-09-27

Branch: `paper-build`, paper directory only.

- Retained the original `paper/paper.pdf`, 48 pages by `pdfinfo`; `pdffonts` shows embedded genuine Times New Roman in that historical artifact. This sandbox lacks the licensed font files and the project's full Pandoc rebuild fails because `xcolor.sty` is unavailable. The original PDF was **not rebuilt** here.
- Added `24_evidence_audit.md`, sourced to committed results and the existing manuscript. Compiled its independent `audit.pdf` (two pages) in the available Times-compatible Nimbus Roman font, not genuine TNR, then concatenated it after the original PDF as `paper_with_audit.pdf` (50 pages by `pdfinfo`). This is a 48+2-page assembled artifact with mixed fonts, not a uniform 50-page TNR paper. Inspected both audit pages; replaced a three-column table that initially collided with a readable numbered claim ledger.
- This audit limits the same-row FoldAmyloid comparison, notes the AntiCP negative, and marks wet-lab and external-validation outcomes pending. No new science result was generated.

## Comparator update, 2026-09-27 morning

- Main now contains committed AntiCP 2.0 comparator identity documentation plus two locked negative iterations. Added `25_comparator_updates.md`, compiled as a two-page `comparator.pdf`, and assembled `paper_with_updates.pdf` from the original 48-page genuine-TNR PDF plus the two-page evidence audit and two-page comparator update. Total 52 pages by `pdfinfo`, but the four appended pages use Nimbus Roman. This is a mixed-font composite, not a uniform-TNR manuscript or a fresh full-paper rebuild.
- Iteration 1 official-test AUROC/MCC 0.7884/0.4192; iteration 2 0.8121/0.4887, both below the published 0.83/0.51 gate. DeepSol and same-partition PASTA 2.0 remain pending. Visually inspected the two new pages; shortened source names and simplified a table after detecting margin overflow.
