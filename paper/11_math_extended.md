# Extended mathematical apparatus

## DeLong variance of an AUROC estimate

Every AUROC in this paper is a Mann-Whitney U statistic. For positives $X_i$ ($i=1..m$) and negatives $Y_j$ ($j=1..n$),

$$\hat{A} = \frac{1}{mn}\sum_{i=1}^{m}\sum_{j=1}^{n} \mathbf{1}[s(X_i) > s(Y_j)] + \tfrac12 \mathbf{1}[s(X_i) = s(Y_j)]$$

Its sampling variance (DeLong et al. 1988) uses the placement values $V_{10}(X_i) = \frac1n \sum_j \psi(X_i, Y_j)$, $V_{01}(Y_j) = \frac1m \sum_i \psi(X_i, Y_j)$ with $\psi(a,b) = \mathbf{1}[a>b] + \frac12\mathbf{1}[a=b]$:

$$\mathrm{Var}(\hat A) = \frac{S_{10}}{m} + \frac{S_{01}}{n}, \qquad S_{10} = \mathrm{Var}_i\big(V_{10}(X_i)\big),\quad S_{01} = \mathrm{Var}_j\big(V_{01}(Y_j)\big)$$

For the pep424 head-to-head ($m=149$, $n=270$, $\hat A = 0.8391$) the DeLong standard error is $\approx 0.020$, so the 95% interval is $[0.800, 0.878]$ - FoldAmyloid's 0.7480 lies more than 4 standard errors below our point estimate. *Derivation sketch.* $\hat A$ is a two-sample U-statistic of degree (1,1); the Hájek projection of a U-statistic onto its first-order terms gives asymptotic normality with variance equal to the variance of the projection, which is exactly the sum of the two placement-value variances above. $\square$

## Metropolis acceptance and why the annealer needs the gates

The screen mutates a sequence $x \to x'$ and accepts with probability

$$p_{\mathrm{acc}} = \min\Big(1,\ \exp\frac{S(x') - S(x)}{T}\Big), \qquad S(x) = w_1 \log P_{acp} + w_2 \log P_{sol} - w_3 \log P_{agg}$$

with $T$ geometrically annealed $T_0 \to T_0 \alpha^k$, $\alpha = 0.95$ per 10 steps. *Claim.* For fixed $T > 0$ the chain's stationary distribution is $\pi_T(x) \propto e^{S(x)/T}$ over the finite space of 20$^L$ sequences. *Proof.* The proposal is symmetric (uniform residue substitution), the chain is irreducible (any sequence reaches any other in $\le L$ mutations) and aperiodic (self-loops via rejection), so detailed balance $\pi(x) p(x \to x') = \pi(x') p(x' \to x)$ - which the Metropolis rule satisfies by construction - implies $\pi_T$ is the unique stationary distribution. $\square$ The corollary is the Section 15 negative: as $T \to 0$ the chain concentrates on $\arg\max S$, and if the learned $P_{acp}$ assigns its maximum to anionic sequences, the chain *correctly* finds the wrong peptides. The gates are not post-hoc filters on a sound optimizer; they are the constraint set that makes the optimization problem the one we mean. Constrained annealing (rejection of gate-violating proposals inside the kernel) is the principled fix and is implemented as the next screen version.

## Additive attention pooling (CNNv2)

For per-position features $h_t \in \mathbb{R}^{c}$ ($t = 1..L$), CNNv2 computes

$$\alpha_t = \frac{\exp(w^\top h_t)}{\sum_{u=1}^{L} \exp(w^\top h_u)}, \qquad z = \sum_{t=1}^{L} \alpha_t h_t$$

*Claim.* $z$ is invariant to zero-padding of the sequence positions when $w^\top h_{\mathrm{pad}} \to -\infty$; with learned scores it is approximately invariant, which is why we mask padded positions before the softmax (implementation detail that changes results: without masking, short peptides in a 60-wide batch get diluted attention, measurably lowering val AUROC in early runs). *Gradient note.* $\partial z / \partial \alpha_t = h_t$ and $\partial \alpha_t / \partial s_t = \alpha_t (1 - \alpha_t)$ with $s_t = w^\top h_t$, so the score path is $O(Lc)$ per example - attention is the cheapest pooling that is content-dependent, which is why it survives our compute budget.

## Early stopping as an $L_2$ shrinkage prior

With quadratic loss $\ell(\theta)$ and Hessian $H$ at the optimum $\theta^*$, stopping gradient descent with step size $\eta$ at iteration $k$ yields parameters $\theta^{(k)} = \big(I - (I - \eta H)^k\big)\,\theta^*$ in the quadratic approximation - a spectral filter that shrinks direction $i$ by factor $1 - (1 - \eta \lambda_i)^k$, exactly the ridge family $\lambda_i / (\lambda_i + \tau)$ with $\tau$ monotone in $1/k$. This is why all our trainers early-stop on validation AUROC rather than training to convergence: the validation carve chooses the effective ridge strength per run, and the honest test number is whatever that choice generalizes to - reported, not tuned.
