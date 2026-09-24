# Data manifest (every entry: source URL verified live on 2026-09-24)

| dataset | role | source | n |
|---|---|---|---|
| AntiCP 2.0 main train/test | ACP benchmark (Agrawal et al. 2021, Brief Bioinform) | https://webs.iiitd.edu.in/raghava/anticp2/download.php (add_data/*_main) | 689+690 train, 171+173 test |
| AntiCP 2.0 alternate train/test | ACP benchmark vs random negatives | same (add_data/*_alternate) | 777+777 train, 193+193 test |
| CancerPPD L-natural | ACP discovery pool (Tyagi et al. 2015, NAR) | https://webs.iiitd.edu.in/raghava/cancerppd/natural_seq_dwn/l_natural.txt | 2,647 canonical ACPs |
| CancerPPD D-natural / non-natural | variant ACPs | same server | 25 / 259+29 |
| eSOL E. coli solubility | solubility benchmark (Niwa et al. 2009/2012) | https://dbarchive.biosciencedbc.jp/data/esol/LATEST/esol.zip | 4,132 measurements |
| UniProt UP000000625 | E. coli K-12 proteome sequences for eSOL mapping | https://rest.uniprot.org/uniprotkb/stream?query=(proteome:UP000000625)&format=fasta | 4,403 |
| AmyloGram/WALTZ-DB full | aggregation benchmark (Kozlowski & Burdukiewicz 2017; Beerten et al. 2015) | https://github.com/michbur/AmyloGramAnalysis (data/amyloid_*_full.fasta) | 421 pos + 1,044 neg |
| AmyloGram/WALTZ-DB benchmark split | external test for aggregation | same (data/amyloid_*_benchmark.fasta) | 269 pos + 746 neg |
| FoldAmyloid predictions on pep424 | published-tool head-to-head | same repo (benchmark/FoldAmyloid_pred.txt) | 419 |
| PASTA 2.0 predictions on pep424 | published-tool head-to-head | same repo (benchmark/pasta2_preds/) | per-protein profiles |
| Guruprasad DIWV scale | instability index | via Biopython Bio.SeqUtils.ProtParamData.DIWV (Guruprasad et al. 1990) | 400 dipeptides |
