# 12. Extended related work

This section surveys, for every external resource this project uses or benchmarks against, what it is, what it contributed, and exactly how it enters our pipeline. It doubles as the audit trail for the tool table (Section 10).

## 12.1 Anticancer-peptide prediction

**AntiCP (Tyagi et al. 2013)** introduced the ACP344/ACP740 benchmarks and showed support vector machines over amino-acid composition separate ACPs from non-ACPs with accuracies around 90% on the easier split. Its binary-motif analysis remains a useful sanity check: ACPs are enriched in lysine/leucine-rich amphipathic stretches. **AntiCP 2.0 (Agrawal et al. 2021)** hardened the task by using antimicrobial-but-not-anticancer peptides as negatives (the "main" set), dropping accuracies into the mid-70s and making the benchmark meaningful; it also shipped the design server we draw the locked splits from. **ACP-MHCNN (2021)** applied a multi-head CNN to the AntiCP 1.0 splits. **mACPpred, iACP, ACPred-FL, ACP-DL** explored feature-fusion and deep variants on the same family of splits. Our reading of this literature: single-split numbers without locked protocols are not comparable, which is why every pepx number states its split, seed, and evaluation discipline.

## 12.2 Protein and peptide solubility

The **eSOL** resource (Niwa et al., PNAS 2009; Nucleic Acids Res 2012) measured the solubility of the *E. coli* K-12 proteome in the PURE cell-free system --- chaperone-free, so the measurement reflects intrinsic sequence properties. **PROSO II** (Smialowski et al. 2012) built a two-layer classifier on a broader but more heterogeneous set. **DeepSol** (Khurana et al. 2018) showed convolutional sequence models beat composition features, reaching ~77% accuracy on an eSOL-derived benchmark. **Protein-Sol** (Hebditch et al. 2017) framed solubility as per-residue regression calibrated to proteome quantiles. **NetSolP** (2022) brought protein language models to the task. We deliberately benchmark on eSOL rather than PROSO II: the measurement is cleaner, the split is ours and locked, and the label (percent soluble) is physical.

## 12.3 Aggregation

**TANGO** (2004) computes statistical-mechanics partition functions over conformational states; **AGGRESCAN** (2007) derives a per-residue aggregation scale from experimental mutagenesis; **PASTA 2.0** (Walsh et al. 2014) models inter-strand pairing energetics; **FoldAmyloid** (2010) uses packing-density expectations; **WALTZ** (2010) combined position-specific scoring on WALTZ-DB hexapeptides; **AmyloGram** (Kozlowski & Burdukiewicz 2017) trained n-gram random forests with the current benchmark record (CV AUROC 0.865 on pep424); **AggreScan3D** extends scoring to structures. The pep424 collection --- 419 experimentally characterized peptides and fragments --- is the common ground on which FoldAmyloid, PASTA 2.0, WALTZ and AmyloGram were all evaluated, which is why we chose it for our head-to-head.

## 12.4 Descriptor scales

Every scale in Section 2 traces to a primary source: Kyte--Doolittle hydropathy (1982, computed from water--vapor transfer free energies); Eisenberg consensus hydrophobicity (1984); Hopp--Woods (1981, antigenicity); Chou--Fasman conformational propensities (1978, statistics over early crystal structures); Zhou aggregation propensity (2004); Guruprasad dipeptide instability (1990, regression against measured half-lives); Boman binding potential (2003); Ikai aliphatic index (1980, thermostability correlate); Bjellqvist pKa set (1993, 2D-gel-calibrated). Using primary scales --- rather than learned embeddings alone --- keeps every model decision physically interpretable and, as Section 4.1 shows, independently auditable.

\newpage

# 13. Architecture details

## 13.1 PeptideCNNv2 forward pass (complete)

Input: index tensor $x \in \{0,\dots,20\}^{B \times L}$ (20 = pad), wide features $w \in \mathbb{R}^{B \times 438}$ standardized by training statistics.

1. Embedding $E \in \mathbb{R}^{21 \times 48}$: $H^0 = E[x]^\top \in \mathbb{R}^{B \times 48 \times L}$.
2. Input projection $H^1 = W_{1\times1} * H^0$, $W_{1\times1} \in \mathbb{R}^{128 \times 48 \times 1}$.
3. Three residual dilated blocks, dilation $d \in \{1,2,4\}$:
$$H^{l+1} = \mathrm{ReLU}\!\left(H^l + W_2^l *_{d} \mathrm{ReLU}(\mathrm{BN}(W_1^l *_{d} H^l))\right)$$
Receptive field after three blocks: $1 + 2\cdot(1+2+4)\cdot(3-1)/2 = 29$ residues, covering typical ACP lengths.
4. Attention pooling (Section 2.9) concatenated with channel-wise max pooling: $z = [c_{attn}; c_{max}] \in \mathbb{R}^{B \times 256}$.
5. Wide fusion: $z' = [z; w] \in \mathbb{R}^{B \times 694}$.
6. Head: $\mathrm{Linear}(694 \to 128)$, ReLU, Dropout(0.25), $\mathrm{Linear}(128 \to 1)$; logit output, BCE loss (Section 2.10).

Parameter count: 196,481. Training on 2 CPU cores: ~60 s for 40 epochs over the 3,897-sequence augmented ACP set.

## 13.2 PeptideGNN graph construction

Nodes: residues. Edges: $\{(i,j) : |i-j| \le 3\}$, motivated by alpha-helical contact order (a helix brings $i, i\pm3, i\pm4$ into contact). Three GCN layers (Section 2.8) with hidden width 64 propagate residue states over this topology; mean pooling over non-pad nodes; fusion with the 18-d descriptor vector; head as above. Parameter count: 16,481 --- deliberately small: the benchmark training sets are 400--3,900 sequences and over-parameterized graph models overfit them immediately (observed: a 4-layer variant lost 3 AUROC points on validation).

## 13.3 TriNet

The PeptideCNNv2 trunk (Section 13.1, 96-channel variant) with three independent two-layer heads. Loss: sum of per-task class-weighted BCEs, each batch drawn from one task round-robin. Descriptor standardization is per-task (Section 9, failure 4). At 25 epochs the shared encoder reaches held-out AUROCs reported in Section 7.1.

## 13.4 Why not a protein language model?

ESM-2/ProtT5 embeddings would likely add 1--3 AUROC on all three benchmarks, but the smallest usable variant (esm2_t6_8M) requires ~90 MB of weights and ~1 GB of RAM at inference for 300-residue inputs --- beyond the 1.9 GB sandbox once training data coexists, and contrary to the program's from-scratch rule for the math. This is recorded as a deliberate scope decision, not an oversight; integrating a distilled embedding is future work.

\newpage

# 14. Reproducibility engineering

## 14.1 Environment

Container: Linux x86_64, 2 CPU, 1.9 GB RAM. Python 3.10.12. Packages: torch 2.14.0+cpu, scikit-learn 1.7.2, numpy 2.2.6, pandas 2.3.3, scipy 1.15.3, matplotlib 3.10.9, biopython 1.88, modlamp 4.3.0, peptides 0.3.2, pyreadr, pytest. Total install: ~250 MB; no GPU used anywhere.

## 14.2 Determinism

All pipelines seed numpy/torch/python-random identically (seed 7; fold offsets 7+f). torch CPU kernels are deterministic at these sizes. The two known non-determinisms: ExtraTrees feature subsampling order under n_jobs=2 (bounded, <0.002 AUROC across reruns) and thread interleaving in the 4-arm ensemble's training order (eliminated by sequential execution).

## 14.3 Data lineage

Every raw file carries its fetch URL and fetch date in data/MANIFEST.md. Every results JSON names its generating script in runs/. Every figure names its source JSON. The paper cites results files, and the repository commit of record is stated in each report --- the chain from claim to bytes is never more than two hops.

\newpage
