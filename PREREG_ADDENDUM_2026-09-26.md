# Pre-registration addendum - 2026-09-26 (rules 1-8 revival)
Locked BEFORE any new outcome is scored. Original registrations unchanged.

## C1 - Benchmark-beat gates (locked)
Verified state (live clone 8a9bdeb1, 29/29 tests green 4:23 PM): cnnv2-hpo test AUC
0.8029 / MCC 0.4638 on anticp2_main (n_test 344) - published AntiCP 2.0 main-dataset
numbers must now be pulled from the paper and the comparison made on the identical
split. Gaps: (A) no same-split comparison vs published AntiCP 2.0 / ACP-MHCNN numbers
yet; (B) solubility arm lacks a same-split DeepSol/Protein-Sol comparison; (C)
aggregation arm lacks a same-split PASTA 2.0 comparison.
Locked protocol: for each of the 3 tasks, rebuild the published benchmark's own split,
score the committed models, and compare against the published numbers on the identical
partition with 10,000-replicate bootstrap CIs. BEAT (per task) = exceed the strongest
published comparator on its primary metric with CI excluding 0. Project benchmark-beat
gate = beat on at least one task with the other two honestly reported; then improve
the remaining arms (feature/architecture iteration, each iteration locked as a new
addendum) until they also beat or a rule-6 pivot redirects them.

## C2 - Tool-gate resolution (locked)
The 40-tool claim was corrected to 24 executed candidates (TOOLS_LEDGER.md). Options,
locked now: (a) execute 16 additional REAL research-tool analyses with committed
evidence (no citation-only counting), or (b) formally restate the gate to the 24-tool
level with her approval via the program lead. No silent redefinition.

## C3 - Discovery arm (unchanged scope)
Tri-objective named candidates re-screened under any improved models; zero-pass =
documented negative + rule-6 pivot.

## Judge rounds
Minimum 10, each producing a concrete novelty improvement (rule 8), verbatim in
JUDGE_ROUNDS.md.
