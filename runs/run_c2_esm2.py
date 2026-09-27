"""C2 ESM-2 embedding baseline (PREREG_C2_ESM2_BASELINE_2026-09-27.md).
Phases: embed (checkpoint) -> score. Frozen esm2_t12_35M_UR50D mean-pooled
embeddings + logistic head per task."""
import sys, json, time
sys.path.insert(0, 'src')
import numpy as np
from pathlib import Path

PART = Path('results/c2_partial')
PART.mkdir(exist_ok=True)

def collect_seqs():
    from pepx.datasets import load_anticp2, load_esol, load_amylogram
    seqs = {}
    tasks = {}
    for split in ('main', 'alternate'):
        recs = load_anticp2(split)
        tasks[f'anticp2_{split}'] = (
            [r.sequence for r in recs if r.meta == 'train'],
            np.array([r.label for r in recs if r.meta == 'train']),
            [r.sequence for r in recs if r.meta == 'test'],
            np.array([r.label for r in recs if r.meta == 'test']))
    esol = load_esol()
    idx = json.load(open('results/esol_split_seed7.json'))
    tr_i, te_i = idx['train_idx'], idx['test_idx']
    tasks['esol'] = ([esol[i].sequence for i in tr_i], np.array([esol[i].label for i in tr_i]),
                     [esol[i].sequence for i in te_i], np.array([esol[i].label for i in te_i]))
    amyl = load_amylogram('full')
    tasks['amylogram_full'] = ([r.sequence for r in amyl], np.array([r.label for r in amyl]), [], np.array([]))
    # pep424 eval sequences
    pep424 = []
    hdr = None
    for line in open('data/raw/aggregation/pep424.fasta'):
        line = line.strip()
        if line.startswith('>'):
            hdr = line[1:]
        elif line and hdr:
            pep424.append(line)
            hdr = None
    lab_by_seq = {}
    for line in open('data/raw/aggregation/pep424_evaluation.txt'):
        p = line.split()
        if len(p) >= 3 and p[-1] in ('+', '-'):
            lab_by_seq[p[1]] = 1 if p[-1] == '+' else 0
    tasks['pep424_raw'] = (pep424, np.array([lab_by_seq.get(s, -1) for s in pep424]), [], np.array([]))
    return tasks

def phase_embed(tasks):
    out_p = PART / 'emb.npz'
    if out_p.exists():
        print('embed done, skip'); return
    import torch, esm
    uniq = sorted({s for t in tasks.values() for s in list(t[0]) + list(t[2])})
    # pep424 hexapeptide windows must be embedded too (window max-pool scoring)
    peps = tasks['pep424_raw'][0]
    wins = {p[j:j+6] for p in peps if set(p) <= set('ACDEFGHIKLMNPQRSTVWY')
            for j in range(len(p) - 5)}
    uniq = sorted(set(uniq) | wins)
    print(f'{len(uniq)} unique sequences', flush=True)
    model, alphabet = esm.pretrained.esm2_t12_35M_UR50D()
    model.eval()
    bc = alphabet.get_batch_converter()
    part_p = PART / 'emb_partial.npz'
    chunks, start = [], 0
    if part_p.exists():
        d = np.load(part_p); chunks = [d['emb']]; start = int(d['done'])
        print('resume at', start, flush=True)
    # token-budget batching (OOM fix): attention peak ~ batch_tokens^2;
    # cap 4096 tokens/batch -> <=~350MB peak for 1k-aa proteins
    # OOM fix v2: attention transient ~ B*H*L_pad^2 (L padded to batch max).
    # Constraint B*Lmax^2 <= 2.5e6 keeps peak <~700MB; length-sorted batches.
    rem = sorted(uniq[start:], key=len)
    batches, cur = [], []
    def lmax(b):
        return min(len(b[-1]), 1022) + 2 if b else 0
    for s in rem:
        B = len(cur) + 1
        L = min(len(s), 1022) + 2
        if cur and (B * L * L > 2_500_000 or B > 128):
            batches.append(cur); cur = []
        cur.append(s)
    if cur:
        batches.append(cur)
    done = start
    ckpt_every = max(1, len(batches) // 40)
    with torch.no_grad():
        for bi, chunk in enumerate(batches):
            _, _, toks = bc([(f'p{i}', p[:1022]) for i, p in enumerate(chunk)])
            rep = model(toks, repr_layers=[12])['representations'][12]
            for r, p in enumerate(chunk):
                L = min(len(p), 1022)
                chunks.append(rep[r, 1:L+1].mean(0).numpy()[None, :])
            done += len(chunk)
            if bi % ckpt_every == 0 or done == len(uniq):
                np.savez_compressed(part_p, emb=np.concatenate(chunks), done=done)
                print(f'embed {done}/{len(uniq)}', flush=True)
                chunks = [np.concatenate(chunks)]
    np.savez_compressed(out_p, seqs=np.array(uniq), emb=np.concatenate(chunks))
    part_p.unlink(missing_ok=True)

def phase_score(tasks):
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score, matthews_corrcoef
    d = np.load(PART / 'emb.npz', allow_pickle=True)
    idx = {p: i for i, p in enumerate(d['seqs'])}
    emb = d['emb']
    def X(seqs):
        return emb[[idx[s] for s in seqs]]
    def boot_ci(fn, y, s, n=10000, seed=23):
        rng = np.random.RandomState(seed); st = []
        y = np.asarray(y); s = np.asarray(s); N = len(y)
        for _ in range(n):
            ii = rng.randint(0, N, N)
            if len(set(y[ii])) < 2:
                continue
            st.append(fn(y[ii], s[ii]))
        return float(np.percentile(st, 2.5)), float(np.percentile(st, 97.5))
    out = {'prereg': 'docs/PREREG_C2_ESM2_BASELINE_2026-09-27.md', 'tasks': {}}
    def decision(v, lo, hi, bar):
        if v > bar and lo > bar: return 'BEAT'
        if lo <= bar <= hi: return 'MATCH'
        return 'MISS'
    for name in ('anticp2_main', 'anticp2_alternate', 'esol'):
        tr_s, tr_y, te_s, te_y = tasks[name]
        clf = LogisticRegression(C=1.0, max_iter=2000).fit(X(tr_s), tr_y)
        s = clf.decision_function(X(te_s))
        p = 1 / (1 + np.exp(-s))
        auc = roc_auc_score(te_y, s)
        mcc = matthews_corrcoef(te_y, (p >= 0.5).astype(int))
        lo_a, hi_a = boot_ci(roc_auc_score, te_y, s)
        lo_m, hi_m = boot_ci(lambda y, q: matthews_corrcoef(y, (q >= 0).astype(int)), te_y, s)
        row = {'auc': float(auc), 'auc_ci95': [lo_a, hi_a], 'mcc': float(mcc),
               'mcc_ci95': [lo_m, hi_m], 'n_test': len(te_y)}
        if name == 'anticp2_main':
            row['comparator'] = {'auc': 0.83, 'mcc': 0.51}
            row['decision_auc'] = decision(auc, lo_a, hi_a, 0.83)
            row['decision_mcc'] = decision(mcc, lo_m, hi_m, 0.51)
        out['tasks'][name] = row
        print(name, row, flush=True)
    # falsifier: shuffled anticp2_main train labels, seed 13
    tr_s, tr_y, te_s, te_y = tasks['anticp2_main']
    rng = np.random.RandomState(13)
    clf = LogisticRegression(C=1.0, max_iter=2000).fit(X(tr_s), rng.permutation(tr_y))
    fals = roc_auc_score(te_y, clf.decision_function(X(te_s)))
    out['falsifier'] = {'seed': 13, 'shuffled_auc': float(fals), 'in_0.40_0.60': bool(0.40 <= fals <= 0.60)}
    print('falsifier', fals, flush=True)
    # pep424: train on amylogram_full hexapeptides, score windows max-pool
    am_s, am_y, _, _ = tasks['amylogram_full']
    clf = LogisticRegression(C=1.0, max_iter=2000).fit(X(am_s), am_y)
    peps, plab, _, _ = tasks['pep424_raw']
    scores, labs = [], []
    win_cache = {}
    for i, pep in enumerate(peps):
        if not set(pep) <= set('ACDEFGHIKLMNPQRSTVWY'):
            continue
        wins = [pep[j:j+6] for j in range(len(pep) - 5)]
        ws = []
        for w in wins:
            if w not in win_cache:
                if w in idx:
                    win_cache[w] = emb[idx[w]]
                else:
                    win_cache[w] = None
            if win_cache[w] is not None:
                ws.append(win_cache[w])
        if not ws:
            continue
        s = clf.decision_function(np.stack(ws)).max()
        scores.append(s)
        lab = plab[i]
        if lab >= 0:
            labs.append(lab)
        else:
            scores.pop()
    auc4 = roc_auc_score(labs, scores)
    lo4, hi4 = boot_ci(roc_auc_score, labs, scores)
    out['tasks']['pep424'] = {'auc': float(auc4), 'auc_ci95': [lo4, hi4], 'n_scored': len(scores),
        'comparators': {'amylogram_published_full424': 0.865, 'foldamyloid_same_subset': 0.5688},
        'decision_vs_same_partition': decision(auc4, lo4, hi4, 0.5688)}
    print('pep424', out['tasks']['pep424'], flush=True)
    json.dump(out, open('results/c2_esm2_baseline.json', 'w'), indent=2)

tasks = collect_seqs()
phase_embed(tasks)
phase_score(tasks)
print('ALL DONE')
