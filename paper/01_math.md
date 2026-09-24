# 2. Mathematical framework

Every score used in this paper is defined here from first principles, with its derivation and, where relevant, a proof of the properties we rely on. Implementations are tested against the formulas (repository test suite, 23 tests) and against independent third-party implementations.

## 2.1 Notation

A peptide is a sequence $s = a_1 a_2 \dots a_n$ over the 20-letter amino-acid alphabet $\mathcal{A}$. $N_g(s)$ counts residues of group $g$ in $s$. Scales are functions $h: \mathcal{A} \to \mathbb{R}$; sequence means are written $\bar{h}(s) = \frac{1}{n}\sum_{i=1}^n h(a_i)$.

## 2.2 Grand average of hydropathicity (F1)

$$\mathrm{GRAVY}(s) = \frac{1}{n}\sum_{i=1}^{n} h_{KD}(a_i)$$

with $h_{KD}$ the Kyte--Doolittle hydropathy scale (1982). Linearity gives $\mathrm{GRAVY}(s_1 s_2) = \frac{n_1 \mathrm{GRAVY}(s_1) + n_2 \mathrm{GRAVY}(s_2)}{n_1 + n_2}$, which we use to verify windowed computations.

## 2.3 Formal charge and the isoelectric point (F2, F3)

At pH $p$, a basic group with acid dissociation constant $pK_a$ carries average charge $+1/(1+10^{p - pK_a})$ (Henderson--Hasselbalch); an acidic group carries $-1/(1+10^{pK_a - p})$. With Bjellqvist et al. (1993) p$K_a$ values, the formal charge is

$$Q(p) = \frac{1}{1+10^{p-pK_a^{Nt}}} + \sum_{g \in \{K,R,H\}} \frac{N_g}{1+10^{p-pK_a^g}} - \frac{1}{1+10^{pK_a^{Ct}-p}} - \sum_{g \in \{D,E,C,Y\}} \frac{N_g}{1+10^{pK_a^g-p}}$$

**Definition (pI).** $pI(s)$ is the unique root of $Q(p) = 0$.

**Proposition (existence, uniqueness, and correctness of bisection).** $Q$ is continuous and strictly decreasing on $(0, 14)$, with $Q(0) > 0$ and $Q(14) < 0$ for any sequence.

*Proof.* Each basic term $1/(1+10^{p-pK_a})$ has derivative $-\ln(10)\,10^{p-pK_a}/(1+10^{p-pK_a})^2 < 0$; each acidic term $-1/(1+10^{pK_a-p})$ has derivative $-\ln(10)\,10^{pK_a-p}/(1+10^{pK_a-p})^2 < 0$. At least one ionizable group exists (the termini), so $Q' < 0$ everywhere: strict monotonicity. As $p \to 0$, basic terms $\to 1$ and acidic terms $\to 0$, so $Q(0) \ge 1$ (N-terminus). As $p \to 14$, basic terms $\to 0$, acidic terms $\to -1$, so $Q(14) \le -1$ (C-terminus). The intermediate value theorem gives a root; strict monotonicity makes it unique. Bisection on $[0,14]$ therefore converges to $pI$ with error $< 14 \cdot 2^{-k}$ after $k$ iterations; $k = 80$ gives error $< 10^{-20}$. $\blacksquare$

## 2.4 Instability index (F4)

Guruprasad, Reddy and Pandit (1990) assign each of the 400 dipeptides an instability weight $\mathrm{DIWV}(xy)$ calibrated on proteins of known half-life:

$$I(s) = \frac{10}{n} \sum_{i=1}^{n-1} \mathrm{DIWV}(a_i a_{i+1})$$

$I < 40$ classifies a protein as stable. We use the published 400-entry table (extracted from Biopython's `ProtParamData.DIWV` and stored in the repository); our implementation reproduces modlAMP's `GlobalDescriptor.instability_index` **exactly** on all test sequences (Section 4, Table 4.1).

## 2.5 Boman index (F5)

$$B(s) = \frac{1}{n}\sum_{i=1}^n b(a_i)$$

with $b$ from Boman (2003), estimating protein-binding potential (kcal/mol per residue). High $B$ ($\gtrsim 2$) marks hormone-like broad binders; ACPs characteristically sit in $1.5 \le B \le 3.2$, a gate we use in the discovery screen. Control: $B(\text{KLAKLAKKLAKLAK}) = 2.561$, matching the published 2.5--2.6 band.

## 2.6 Aliphatic index (F6)

$$\mathrm{AI}(s) = X_A + 2.9\,X_I + 3.9\,(X_L + X_V)$$

with $X_g$ the mole percent of residue $g$ (Ikai 1980). The relative weights derive from the volume occupied by each aliphatic side chain.

## 2.7 Amino-acid and dipeptide composition (F7)

$$\mathrm{AAC}_g(s) = \frac{N_g(s)}{n}, \qquad \mathrm{DPC}_{gh}(s) = \frac{|\{i : a_i = g, a_{i+1} = h\}|}{n-1}$$

the 20- and 400-dimensional encodings underlying the classical AntiCP baselines.

## 2.8 Graph convolution on the sequence chain (F8)

Residues are nodes; edges join residues within sequence distance $w$ (we use $w = 3$): $A_{ij} = \mathbf{1}[|i-j| \le w]$. With $\tilde{A} = A + I$ and $\tilde{D}_{ii} = \sum_j \tilde{A}_{ij}$,

$$H^{(l+1)} = \sigma\!\left(\tilde{D}^{-1/2} \tilde{A} \tilde{D}^{-1/2} H^{(l)} W^{(l)}\right)$$

The symmetric normalization keeps eigenvalues of the propagation operator in $[-1, 1]$: $\tilde{D}^{-1/2}\tilde{A}\tilde{D}^{-1/2} = I - \tilde{L}_{sym}$ with $\tilde{L}_{sym}$ the normalized Laplacian of $\tilde{A}$, whose spectrum lies in $[0,2]$. This bound is what stabilizes repeated application (Kipf & Welling 2017). Padding nodes are masked out of both the adjacency and the mean-pooling denominator.

## 2.9 Additive attention pooling (F9)

$$\alpha_i = \frac{\exp(w^\top h_i)}{\sum_{j \in R} \exp(w^\top h_j)}, \qquad c = \sum_{i \in R} \alpha_i h_i$$

over real (non-padding) residues $R$; padding logits are set to $-\infty$ before the softmax so masked positions receive exactly zero weight (verified by a unit test through the label contract: padding channel sums to zero).

## 2.10 Class-weighted binary cross-entropy (F10)

$$\mathcal{L} = -\frac{1}{N}\sum_{k=1}^N \left[ w_p\, y_k \log \sigma(z_k) + (1 - y_k) \log (1 - \sigma(z_k)) \right], \quad w_p = \frac{N_{neg}}{N_{pos}}$$

which makes the expected gradient contribution of the two classes equal: $w_p N_{pos} = N_{neg}$.

## 2.11 Tri-objective screen score and annealing (F11)

$$\mathrm{S}(s) = P_{acp}(s) \cdot P_{sol}(s) \cdot \left(1 - P_{agg}(s)\right)$$

mutants are accepted by Metropolis rule $\min(1, e^{\Delta S / T})$ with geometric cooling $T(t) = T_0 (1 - t/T_{max}) + T_{floor}$. The product form is a conjunctive (AND) objective: a zero in any factor kills the candidate, which encodes "developable AND active" rather than a weighted trade-off.

## 2.12 Evaluation metrics (F12)

$$\mathrm{MCC} = \frac{TP \cdot TN - FP \cdot FN}{\sqrt{(TP+FP)(TP+FN)(TN+FP)(TN+FN)}}$$

$$\mathrm{AUROC} = P\left(s(x^+) > s(x^-)\right) = \frac{1}{N_+ N_-} \sum_{i,j} \mathbf{1}[s(x_i^+) > s(x_j^-)] + \tfrac{1}{2}\mathbf{1}[s(x_i^+) = s(x_j^-)]$$

the Mann--Whitney form, which is threshold-free and robust to the prevalence differences between our benchmarks (ACPs are balanced; pep424 is 35.6% positive; eSOL-at-30% is 62.3% positive).

\newpage
