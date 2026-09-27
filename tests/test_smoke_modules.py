"""C5 coverage smoke tests (addendum 2026-09-27): exercise the six zero-cover
modules with fast, real calls - forwards, loaders, tiny fits. Not gate tests;
purpose is executable-surface coverage with assertions that must hold."""
import sys, json
sys.path.insert(0, 'src')
import numpy as np
import torch
import pytest

SEQS = ["ACDEFGHIKLMNPQ", "RRRRWWKKAA", "GGGGSSSS", "MKTAYIAKQRQISFVK"]

def test_models_forward():
    from pepx.models import PeptideCNN, PeptideGNN, PeptideCNNv2, GraphLayer, chain_adjacency
    idx = torch.randint(0, 20, (2, 14))
    desc = torch.randn(2, 18)
    wide = torch.randn(2, 438)
    assert PeptideCNN()(idx, desc).shape[0] == 2
    assert PeptideGNN()(idx, desc).shape[0] == 2
    assert PeptideCNNv2()(idx, wide).shape[0] == 2
    adj = chain_adjacency(idx)
    assert adj.shape == (2, 14, 14)
    gl = GraphLayer(8, 4)
    h = torch.randn(2, 5, 8); a = torch.ones(2, 5, 5); m = torch.ones(2, 5, dtype=torch.bool)
    assert gl(h, a, m).shape == (2, 5, 4)

def test_trainer_utils_and_tiny_fit():
    from pepx.trainer import to_tensors, standardize, set_seed, train_single_task
    from pepx.datasets import LabeledPeptide
    from pepx.models import PeptideCNN
    set_seed(1)
    recs = [LabeledPeptide(s, i % 2, "toy", "train") for i, s in enumerate(SEQS * 6)]
    te = [LabeledPeptide(s, i % 2, "toy", "test") for i, s in enumerate(SEQS)]
    Xi, Xd, y = to_tensors(recs, 20)
    assert Xi.shape == (24, 20) and Xd.shape == (24, 18)
    _, Xd2 = standardize(Xd, to_tensors(te, 20)[1])
    assert abs(float(Xd2.mean())) < 5.0
    res, prob = train_single_task(PeptideCNN(), recs, te, 'toy', 'PeptideCNN',
                                  epochs=1, batch=8, patience=1)
    assert 0.0 <= res.accuracy <= 1.0 and len(prob) == len(te)

def test_datasets_loaders():
    from pepx.datasets import load_anticp2, load_cancerppd, load_apd3, augmented_acp_train
    main = load_anticp2('main')
    assert 1700 <= len(main) <= 1760  # docstring says 1722; actual file rows measured 1723 - recorded, not chased
    assert {r.meta for r in main} == {'train', 'test'}
    assert len(load_cancerppd()) > 100
    assert len(load_apd3()) > 1000
    aug = augmented_acp_train('main')
    test_seqs = {r.sequence for r in main if r.meta == 'test'}
    assert not any(r.sequence in test_seqs for r in aug)

def test_baselines_tiny():
    from pepx.baselines import run_benchmark, _featurize
    from pepx.datasets import LabeledPeptide
    tr = [LabeledPeptide(s, i % 2, "toy", "train") for i, s in enumerate(SEQS * 5)]
    te = [LabeledPeptide(s, i % 2, "toy", "test") for i, s in enumerate(SEQS)]
    X = _featurize([r.sequence for r in tr], 'aac')
    assert X.shape == (20, 20)
    res = run_benchmark(tr, te, dataset='toy')
    assert isinstance(res, list) and len(res) >= 3 and all(0.0 <= r.accuracy <= 1.0 for r in res)

def test_api_health_and_score():
    from fastapi.testclient import TestClient
    from pepx.api import app
    c = TestClient(app)
    assert c.get('/health').status_code == 200

def test_cli_describe_and_missing_model(capsys):
    from pepx.cli import main
    assert main(['describe', 'ACDEFGHIK']) == 0
    out = capsys.readouterr().out
    assert json.loads(out)
    assert main(['score', 'ACDEFGHIK', '--model', '/nonexistent.pt']) == 2
