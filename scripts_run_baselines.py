"""Run classical baselines on the AntiCP 2.0 main/alternate locked splits."""
import sys, json
sys.path.insert(0, 'src')
from pepx.datasets import load_anticp2
from pepx.baselines import run_benchmark

results = []
for split in ("main", "alternate"):
    recs = load_anticp2(split)
    train = [r for r in recs if r.meta == "train"]
    test = [r for r in recs if r.meta == "test"]
    results += run_benchmark(train, test, f"anticp2_{split}")

for r in results:
    print(f"{r.dataset:20s} {r.model:8s} acc={r.accuracy:.4f} sn={r.sensitivity:.4f} "
          f"sp={r.specificity:.4f} mcc={r.mcc:.4f} auc={r.auc:.4f}")
with open("results/baselines_anticp2.json", "w") as fh:
    json.dump([r.__dict__ for r in results], fh, indent=2)
