"""Classical published-baseline reproductions: AAC/DPC + SVM / RF / ET,
the exact baseline stack of AntiCP/AntiCP 2.0. These are the floor every
deep model must beat."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, List, Tuple

import numpy as np
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.metrics import (accuracy_score, matthews_corrcoef, roc_auc_score,
                             confusion_matrix)
from sklearn.svm import SVC

from .datasets import LabeledPeptide
from .encoders import aac, dpc, descriptor_array


@dataclass
class BenchmarkResult:
    model: str
    encoding: str
    dataset: str
    accuracy: float
    sensitivity: float
    specificity: float
    mcc: float
    auc: float
    n_train: int
    n_test: int


def _featurize(seqs: List[str], enc: str) -> np.ndarray:
    f: Callable[[str], np.ndarray] = {"aac": aac, "dpc": dpc, "desc": descriptor_array}[enc]
    return np.stack([f(s) for s in seqs])


def run_benchmark(train: List[LabeledPeptide], test: List[LabeledPeptide],
                  dataset: str) -> List[BenchmarkResult]:
    """Train the classical grid on train, evaluate once on the locked test split."""
    Xtr = {"aac": _featurize([r.sequence for r in train], "aac"),
           "dpc": _featurize([r.sequence for r in train], "dpc"),
           "desc": _featurize([r.sequence for r in train], "desc")}
    Xte = {"aac": _featurize([r.sequence for r in test], "aac"),
           "dpc": _featurize([r.sequence for r in test], "dpc"),
           "desc": _featurize([r.sequence for r in test], "desc")}
    ytr = np.array([r.label for r in train])
    yte = np.array([r.label for r in test])
    models: Dict[str, Tuple[str, object]] = {
        "SVM-AAC": ("aac", SVC(kernel="rbf", C=10, gamma=0.01, probability=True, random_state=7)),
        "SVM-DPC": ("dpc", SVC(kernel="rbf", C=10, gamma=0.001, probability=True, random_state=7)),
        "RF-DPC": ("dpc", RandomForestClassifier(n_estimators=300, random_state=7, n_jobs=2)),
        "ET-DPC": ("dpc", ExtraTreesClassifier(n_estimators=300, random_state=7, n_jobs=2)),
        "RF-DESC": ("desc", RandomForestClassifier(n_estimators=300, random_state=7, n_jobs=2)),
    }
    out: List[BenchmarkResult] = []
    for name, (enc, clf) in models.items():
        clf.fit(Xtr[enc], ytr)
        pred = clf.predict(Xte[enc])
        prob = clf.predict_proba(Xte[enc])[:, 1]
        tn, fp, fn, tp = confusion_matrix(yte, pred, labels=[0, 1]).ravel()
        out.append(BenchmarkResult(
            model=name, encoding=enc, dataset=dataset,
            accuracy=accuracy_score(yte, pred),
            sensitivity=tp / max(tp + fn, 1),
            specificity=tn / max(tn + fp, 1),
            mcc=matthews_corrcoef(yte, pred),
            auc=roc_auc_score(yte, prob),
            n_train=len(train), n_test=len(test)))
    return out
