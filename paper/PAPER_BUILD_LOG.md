# Paper-build evidence log, 2026-09-27

Branch: `paper-build`, paper directory only.

- Retained the original `paper/paper.pdf`, 48 pages by `pdfinfo`; `pdffonts` shows embedded genuine Times New Roman in that historical artifact. This sandbox lacks the licensed font files and the project's full Pandoc rebuild fails because `xcolor.sty` is unavailable. The original PDF was **not rebuilt** here.
- Added `24_evidence_audit.md`, sourced to committed results and the existing manuscript. Compiled its independent `audit.pdf` (two pages) in the available Times-compatible Nimbus Roman font, not genuine TNR, then concatenated it after the original PDF as `paper_with_audit.pdf` (50 pages by `pdfinfo`). This is a 48+2-page assembled artifact with mixed fonts, not a uniform 50-page TNR paper. Inspected both audit pages; replaced a three-column table that initially collided with a readable numbered claim ledger.
- This audit limits the same-row FoldAmyloid comparison, notes the AntiCP negative, and marks wet-lab and external-validation outcomes pending. No new science result was generated.

## Comparator update, 2026-09-27 morning

- Main now contains committed AntiCP 2.0 comparator identity documentation plus two locked negative iterations. Added `25_comparator_updates.md`, compiled as a two-page `comparator.pdf`, and assembled `paper_with_updates.pdf` from the original 48-page genuine-TNR PDF plus the two-page evidence audit and two-page comparator update. Total 52 pages by `pdfinfo`, but the four appended pages use Nimbus Roman. This is a mixed-font composite, not a uniform-TNR manuscript or a fresh full-paper rebuild.
- Iteration 1 official-test AUROC/MCC 0.7884/0.4192; iteration 2 0.8121/0.4887, both below the published 0.83/0.51 gate. DeepSol and same-partition PASTA 2.0 remain pending. Visually inspected the two new pages; shortened source names and simplified a table after detecting margin overflow.

## Judge requirement amended, 2026-09-27 10:00 IST

The owner said "NOT 10 ROUNDS OOF CHATGPT CHECK JUST ONE WHICH I PROVIDE OK?" (authenticated WhatsApp message `wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMDJCMTZGRTVEMkQwMTFBQzc4MQA=`, 10:00:07 IST). For this project, the paper branch therefore marks the counted judge gate **0 of 1, PENDING her personally provided verdict**. A round initiated by agents, even through her ChatGPT account, remains historical or supplementary and does not meet the gate. Historical ten-round language in the science ledgers is not erased by this note. A project-specific user-pasted verdict must be traced and evaluated before completion is recorded. Supplementary Gemini/LLM consults do not count. No scientific result, page count or font gate changes here.

## Authorship-attribution cleanup, 2026-09-27 11:14 IST

The owner requested removal of the assistant's attribution from the papers (WhatsApp `wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMEY5MzY4M0Q4OUYwNjg4ODZDNwA=`). Removed agent/program-style byline and credit text from the editable paper source and PDF display, without substituting an author. Udita's own byline in 09b was preserved, with only the Instinct pipeline parenthetical removed. Manuscript PDF author metadata is empty. Literature references to other studies' authors and technical uses of "author numbering" are not authorship credits for this paper.

Licensed-TNR rebuild still unavailable here: the three existing paper PDFs (48-page source PDF, 50-page audit composite, 52-page updates composite) were redacted on the first page to remove the generated pipeline byline. This is a PDF edit, not a fresh typeset rebuild. First-page pixels and text extraction were checked; licensed TNR remains embedded on the original 48 pages, while appended pages in the composites remain Nimbus Roman.
