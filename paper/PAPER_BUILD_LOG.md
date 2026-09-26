# Paper-build evidence log, 2026-09-27

Branch: `paper-build`, paper directory only.

- Retained the original `paper/paper.pdf`, 48 pages by `pdfinfo`; `pdffonts` shows embedded genuine Times New Roman in that historical artifact. This sandbox lacks the licensed font files and the project's full Pandoc rebuild fails because `xcolor.sty` is unavailable. The original PDF was **not rebuilt** here.
- Added `24_evidence_audit.md`, sourced to committed results and the existing manuscript. Compiled its independent `audit.pdf` (two pages) in the available Times-compatible Nimbus Roman font, not genuine TNR, then concatenated it after the original PDF as `paper_with_audit.pdf` (50 pages by `pdfinfo`). This is a 48+2-page assembled artifact with mixed fonts, not a uniform 50-page TNR paper. Inspected both audit pages; replaced a three-column table that initially collided with a readable numbered claim ledger.
- This audit limits the same-row FoldAmyloid comparison, notes the AntiCP negative, and marks wet-lab and external-validation outcomes pending. No new science result was generated.
