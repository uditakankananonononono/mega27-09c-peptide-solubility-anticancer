"""eSOL solubility benchmark: locked stratified split (seed 7), classical +
deep models. Reference: DeepSol (Khurana et al. 2018) ~77% acc on its
eSOL-derived benchmark."""
import sys, json, time
sys.path.insert(0, 'src')
import numpy as np
from sklearn.model_selection import train_test_split
from pepx.datasets import load_esol
from pepx.baselines import run_benchmark
from pepx.models import PeptideCNNv2, PeptideGNN
from pepx.trainer import train_single_task

t0 = time.time()
recs = load_esol()
y = np.array([r.label for r in recs])
idx_tr, idx_te = train_test_split(np.arange(len(recs)), test_size=0.15,
                                  stratify=y, random_state=7)
train = [recs[i] for i in idx_tr]
test = [recs[i] for i in idx_te]
json.dump({'train_idx': idx_tr.tolist(), 'test_idx': idx_te.tolist()},
          open('results/esol_split_seed7.json', 'w'))
print(f"train={len(train)} test={len(test)}", flush=True)

for r in run_benchmark(train, test, 'esol'):
    print(f"{r.model:8s} acc={r.accuracy:.4f} mcc={r.mcc:.4f} auc={r.auc:.4f}", flush=True)

res, _ = train_single_task(PeptideCNNv2(), train, test, 'esol', 'CNNv2',
                           max_len=300, epochs=12, wide=True)
print(res, flush=True)
res2, _ = train_single_task(PeptideGNN(), train, test, 'esol', 'GNN',
                            max_len=300, epochs=12)
print(res2, flush=True)
print('elapsed', round(time.time() - t0, 1), 's', flush=True)
