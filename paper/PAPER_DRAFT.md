# Tri-Objective Peptide Intelligence: Benchmark-Headed Models for Anticancer Activity, Solubility and Aggregation, and a Multi-Task Discovery Screen

MEGA-PROGRAM-27, item 9 part 3. Draft skeleton - Times New Roman target, 50 pp expanded at ship time.

## 1. Abstract
## 2. Introduction (ACP therapeutics; developability triad: activity, solubility, aggregation)
## 3. Mathematical framework (10+ formulas, full derivations in Appendix A)
  F1. Grand average of hydropathicity: GRAVY = (1/n) Σ_i h(a_i)  [Kyte-Doolittle 1982]
  F2. Net charge at pH p: Q(p) = Σ_g N_g/(1+10^(p-pKa_g)) - Σ_h N_h/(1+10^(pKa_h-p)) + termini  [Henderson-Hasselbalch]
  F3. Isoelectric point: pI = arg_p { Q(p) = 0 } by bisection (proof of bracketing: Q monotone decreasing, Q(0)>0>Q(14))
  F4. Guruprasad instability: I = (10/n) Σ_{i=1}^{n-1} DIWV(x_i x_{i+1})  [Protein Eng 1990]
  F5. Boman index: B = (1/n) Σ_i b(a_i)  [AAC 2003]
  F6. Aliphatic index: AI = X(A) + 2.9 X(I) + 3.9 (X(L)+X(V))  [Ikai 1980]
  F7. GCN layer: H^(l+1) = σ( D^{-1/2} (A+I) D^{-1/2} H^(l) W^(l) )  [Kipf & Welling 2017] + derivation of normalization
  F8. Additive attention pooling: α = softmax(W_a h_i); c = Σ α_i h_i
  F9. Class-weighted BCE: L = -[ w_p y log σ(z) + (1-y) log(1-σ(z)) ], w_p = N_neg/N_pos
  F10. Tri-objective screen score: S(s) = P_acp(s) · P_sol(s) · (1 - P_agg(s)), with simulated-annealing acceptance min(1, exp(ΔS/T))
  F11. MCC = (TP·TN - FP·FN) / sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN))
  F12. Dilated residual conv block: y = relu(x + W2 * relu(BN(W1 *_d x)))
## 4. Data (manifest table; decontamination protocol; locked splits)
## 5. Models (classical grid; PeptideCNN; PeptideGNN; CNNv2; TriNet)
## 6. Benchmark results (head-to-head tables vs AntiCP 2.0, DeepSol, FoldAmyloid, PASTA 2.0)
## 7. Discovery screen (named candidates, gates, falsifiability protocol)
## 8. External tools used (40-row table: tool, version, what it was used for here)
## 9. Honest negatives and failure modes
## 10. Conclusion
## Appendix A: derivations and proofs
## Appendix B: full dataset manifest (120+ accession-level sources)
## Appendix C: reproduction (seeds, commands, environment)
