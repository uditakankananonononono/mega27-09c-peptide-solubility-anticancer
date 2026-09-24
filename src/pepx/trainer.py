"""Training harness: fixed seeds, stratified validation carve, early stop on
validation AUC, single final evaluation on the locked test split."""
from __future__ import annotations

import random
from dataclasses import dataclass, asdict
from typing import Dict, List, Tuple

import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import (accuracy_score, confusion_matrix,
                             matthews_corrcoef, roc_auc_score)
from torch.utils.data import DataLoader, TensorDataset

from .alphabet import AA_STANDARD
from .datasets import LabeledPeptide
from .encoders import descriptor_array

PAD = len(AA_STANDARD)
_AA = {a: i for i, a in enumerate(AA_STANDARD)}


def set_seed(seed: int = 7) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def to_tensors(recs: List[LabeledPeptide], max_len: int, wide: bool = False) -> Tuple[torch.Tensor, ...]:
    n = len(recs)
    idx = torch.full((n, max_len), PAD, dtype=torch.long)
    desc = torch.zeros((n, 438 if wide else 18), dtype=torch.float32)
    y = torch.zeros(n, dtype=torch.float32)
    from .encoders import aac, dpc
    for i, r in enumerate(recs):
        ids = [_AA.get(a, PAD) for a in r.sequence[:max_len]]
        idx[i, : len(ids)] = torch.tensor(ids)
        if wide:
            desc[i] = torch.cat([torch.from_numpy(descriptor_array(r.sequence)),
                                 torch.from_numpy(aac(r.sequence)),
                                 torch.from_numpy(dpc(r.sequence))])
        else:
            desc[i] = torch.from_numpy(descriptor_array(r.sequence))
        y[i] = float(r.label)
    # standardize descriptors on this split (train stats applied by caller via fit=False trick:
    # we standardize per-tensor here; splits share stats via standardize_with)
    return idx, desc, y


def standardize(train_desc: torch.Tensor, *others: torch.Tensor) -> Tuple[torch.Tensor, ...]:
    mu = train_desc.mean(0, keepdim=True)
    sd = train_desc.std(0, keepdim=True).clamp(min=1e-6)
    out = [(train_desc - mu) / sd] + [(o - mu) / sd for o in others]
    return tuple(out)


@dataclass
class EvalResult:
    model: str
    dataset: str
    accuracy: float
    sensitivity: float
    specificity: float
    mcc: float
    auc: float
    epochs_run: int
    best_val_auc: float


def _metrics(y_true: np.ndarray, prob: np.ndarray) -> Dict[str, float]:
    pred = (prob >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, pred, labels=[0, 1]).ravel()
    return dict(accuracy=float(accuracy_score(y_true, pred)),
                sensitivity=float(tp / max(tp + fn, 1)),
                specificity=float(tn / max(tn + fp, 1)),
                mcc=float(matthews_corrcoef(y_true, pred)),
                auc=float(roc_auc_score(y_true, prob)))


def train_single_task(model: nn.Module, train: List[LabeledPeptide],
                      test: List[LabeledPeptide], dataset: str, model_name: str,
                      max_len: int = 60, epochs: int = 40, batch: int = 64,
                      lr: float = 1e-3, val_frac: float = 0.15, wide: bool = False,
                      patience: int = 8, seed: int = 7,
                      task: str | None = None) -> Tuple[EvalResult, np.ndarray]:
    """Full disciplined loop. Returns (metrics, test_probabilities)."""
    set_seed(seed)
    y_all = np.array([r.label for r in train])
    rng = np.random.RandomState(seed)
    pos_idx = np.where(y_all == 1)[0]; neg_idx = np.where(y_all == 0)[0]
    rng.shuffle(pos_idx); rng.shuffle(neg_idx)
    nvp = max(1, int(len(pos_idx) * val_frac)); nvn = max(1, int(len(neg_idx) * val_frac))
    val_set = set(pos_idx[:nvp].tolist()) | set(neg_idx[:nvn].tolist())
    tr = [r for i, r in enumerate(train) if i not in val_set]
    va = [r for i, r in enumerate(train) if i in val_set]

    Xtr_i, Xtr_d, ytr = to_tensors(tr, max_len, wide)
    Xva_i, Xva_d, yva = to_tensors(va, max_len, wide)
    Xte_i, Xte_d, yte = to_tensors(test, max_len, wide)
    Xtr_d, Xva_d, Xte_d = standardize(Xtr_d, Xva_d, Xte_d)

    loader = DataLoader(TensorDataset(Xtr_i, Xtr_d, ytr), batch_size=batch, shuffle=True)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    lossf = nn.BCEWithLogitsLoss(pos_weight=torch.tensor(
        max((ytr == 0).sum().item(), 1.0) / max((ytr == 1).sum().item(), 1.0)))

    def forward(i, d):
        return model(i, d, task) if task else model(i, d)

    best_auc, best_state, wait, run = -1.0, None, 0, 0
    for ep in range(1, epochs + 1):
        run = ep
        model.train()
        for i, d, y in loader:
            opt.zero_grad()
            loss = lossf(forward(i, d), y)
            loss.backward()
            opt.step()
        model.eval()
        with torch.no_grad():
            pv = torch.sigmoid(forward(Xva_i, Xva_d)).numpy()
        va_auc = roc_auc_score(yva.numpy(), pv)
        if va_auc > best_auc:
            best_auc, best_state, wait = va_auc, {k: v.clone() for k, v in model.state_dict().items()}, 0
        else:
            wait += 1
            if wait >= patience:
                break
    if best_state:
        model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        pt = torch.sigmoid(forward(Xte_i, Xte_d)).numpy()
    m = _metrics(yte.numpy().astype(int), pt)
    return (EvalResult(model=model_name, dataset=dataset, epochs_run=run,
                       best_val_auc=float(best_auc), **m), pt)
