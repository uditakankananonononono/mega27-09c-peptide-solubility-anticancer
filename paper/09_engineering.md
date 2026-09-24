# Engineering for use: acceleration, portability, serving

A benchmark result is only as useful as the artifact that carries it. This section documents the three engineering steps that turn the research code into a deployable tool, plus one integrity incident we consider part of the result.

## Kernel acceleration with numba

The discovery screen's inner loop - score, mutate, accept/reject over millions of steps - is pure Python in its reference form. We JIT-compiled the accept/reject kernel with **numba** (`@njit`, cache enabled) and benchmarked 2,000,000 steps: **80x speedup** over the interpreted loop (`results/tool_engineering.json`; wall times recorded per run, 80.4x on the committed run). The annealer's semantics are unchanged; the kernel is bit-identical in accept/reject decisions for a fixed seed because no floating-point reordering occurs inside the loop body.

## Model portability with ONNX

The shipped TriNet is exported to ONNX (opset 17, `results/trinet.onnx`, 434 KB) via an all-heads wrapper (the research forward takes a task string, which is not traceable; the wrapper computes all three heads in one graph). Parity against the PyTorch reference on the gated candidate, with descriptors standardized by the shipped train statistics exactly as the CLI does: **max absolute logit difference 4.8e-7** (`results/tool_engineering.json`). The exported graph reproduces the paper's candidate numbers (P_acp 0.9033, P_sol 0.8667, P_agg 0.0550), so downstream consumers do not need a PyTorch install.

## Serving with FastAPI

`pepx.api` exposes the tool over HTTP: `GET /health` (liveness + model version) and `POST /score` (batch up to 512 sequences, per-sequence validation with explicit rejection of non-standard residues). The smoke test (`results/api_smoke.json`) runs the FastAPI TestClient against the app and asserts exact agreement with the CLI on the gated candidate (0.9033) and the KLAKLAK control (0.9899). Serving adds no new model risk: it is the same weights, same normalization, same code path as the CLI, with pydantic request validation at the boundary.

## Class-balance check with imbalanced-learn

Because the AntiCP 2.0 main split is already balanced (689/689 train), applying SMOTE should be neutral - and it is: logistic baseline 0.6339 AUROC with SMOTE vs 0.6337 without (`results/tool_engineering.json`). We record this as a negative control that validates the split construction rather than as a modeling result.

## Integrity incident: the shipped weights were almost silently degraded

During the scaled screen (Section 15), `run_discovery_scaled.py` wrote its 25-epoch weights to `results/trinet.pt` - the same path the CLI, API and ONNX export ship. The overwrite was caught before being committed, but only because the ONNX parity demo returned saturated probabilities: the 25-epoch head, which Section 15 shows drifted toward anionic sequences, scored the gated cationic candidate at P_acp 0.315 instead of 0.904. We restored the v0.1 weights from git (the drifted weights were never committed), preserved them as `results/trinet_25ep.pt` for the negative result, and changed the scaled script to write `_25ep`-suffixed paths. Three lessons, all now enforced: (i) shipped artifacts and experiment artifacts must never share a path; (ii) a tool's smoke test is a model-regression test - the parity demo caught what the training logs did not flag; (iii) "the CLI reproduces the paper's numbers" is a falsifiable claim we re-check after every model change, and it is the reason the discovery claim and the shipped tool cannot silently diverge.
