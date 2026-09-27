"""Gate-constrained annealing screen (the Section 15 fix, implemented):
loads the SHIPPED v0.1 TriNet and mutates only within the descriptor-gate
region - gate-violating proposals are rejected before scoring (the
constrained-kernel variant of the Metropolis chain)."""
import sys, json, time
sys.path.insert(0, 'src')
import numpy as np
import torch
from pepx.datasets import load_anticp2, load_cancerppd, load_apd3
from pepx.models import TriNet
from pepx.trainer import to_tensors, set_seed
from pepx.descriptors import descriptor_vector
from pepx.alphabet import AA_STANDARD

t0 = time.time()
set_seed(7)
known = set(load_cancerppd()) | set(load_apd3())
known |= {r.sequence for r in load_anticp2('main')} | {r.sequence for r in load_anticp2('alternate')}

model = TriNet()
model.load_state_dict(torch.load('results/trinet.pt', map_location='cpu'))
model.eval()
z = np.load('results/trinet_norm.npz'); MU, SD = z['mu'], z['sd']

def tri_score(seq):
    i, d, _ = to_tensors([type('R', (), {'sequence': seq, 'label': 0, 'meta': ''})()], 60)
    d = (d - torch.from_numpy(MU).float()) / torch.from_numpy(SD).float()
    with torch.no_grad():
        s = model.score_all(i, d)
    pa, ps, pg = float(s['acp'][0]), float(s['sol'][0]), float(s['agg'][0])
    return pa * ps * (1 - pg), pa, ps, pg

def gates_ok(seq):
    dv = descriptor_vector(seq)
    return (2 <= dv['net_charge'] <= 9 and dv['gravy'] < 0 and 1.5 <= dv['boman'] <= 3.2
            and seq.count('C') % 2 == 0)  # C3: even-Cys disulfide-aware gate

rng = np.random.RandomState(109)
cands = {}
seeds = [s for s in load_cancerppd() if 13 <= len(s) <= 16][:60]
n_prop, n_rej = 0, 0
for seed in seeds:
    cur = list(seed)
    if not gates_ok(''.join(cur)): continue
    cur_s, *_ = tri_score(''.join(cur))
    for step in range(60):
        pos = rng.randint(len(cur))
        mut = cur.copy()
        mut[pos] = AA_STANDARD[rng.randint(20)]
        mseq = ''.join(mut)
        n_prop += 1
        if not gates_ok(mseq):       # constrained kernel: reject before scoring
            n_rej += 1
            continue
        s_mut, *_ = tri_score(mseq)
        T = 0.05 * (1 - step / 90) + 0.005
        if s_mut > cur_s or rng.rand() < np.exp((s_mut - cur_s) / T):
            cur, cur_s = mut, s_mut
    seq = ''.join(cur)
    if seq not in known:
        sc, pa, ps, pg = tri_score(seq)
        dv = descriptor_vector(seq)
        cands[seq] = dict(score=sc, p_acp=pa, p_sol=ps, p_agg=pg,
                          charge=dv['net_charge'], gravy=dv['gravy'], boman=dv['boman'])

out = []
for seq, m in sorted(cands.items(), key=lambda kv: -kv[1]['score']):
    gates = (m['p_acp'] >= 0.9 and m['p_sol'] >= 0.7 and m['p_agg'] <= 0.2
             and 2 <= m['charge'] <= 9 and m['gravy'] < 0 and 1.5 <= m['boman'] <= 3.2)
    out.append({'sequence': seq, **m, 'passes_all_gates': bool(gates)})
# ---- preregistered controls (C7): scramble (seed 103) + random panel (seed 107)
def gates_tuple(m):
    return (m['p_acp'] >= 0.9 and m['p_sol'] >= 0.7 and m['p_agg'] <= 0.2
            and 2 <= m['charge'] <= 9 and m['gravy'] < 0 and 1.5 <= m['boman'] <= 3.2)

def score_seq(seq):
    sc, pa, ps, pg = tri_score(seq)
    dv = descriptor_vector(seq)
    return dict(score=sc, p_acp=pa, p_sol=ps, p_agg=pg,
                charge=dv['net_charge'], gravy=dv['gravy'], boman=dv['boman'])

passing = [o for o in out if o['passes_all_gates']]
rng_s = np.random.RandomState(111)
scrambles = []
for o in passing:
    ch = list(o['sequence']); rng_s.shuffle(ch); sq = ''.join(ch)
    m = score_seq(sq)
    scrambles.append({'sequence': sq, 'parent': o['sequence'], **m,
                      'passes_all_gates': bool(gates_tuple(m) and gates_ok(sq))})
rng_r = np.random.RandomState(113)
panel = []
for _ in range(100):
    L = int(rng_r.randint(13, 17))
    sq = ''.join(AA_STANDARD[rng_r.randint(20)] for _ in range(L))
    if not gates_ok(sq):
        panel.append({'sequence': sq, 'in_kernel': False, 'passes_all_gates': False})
        continue
    m = score_seq(sq)
    panel.append({'sequence': sq, 'in_kernel': True, **m,
                  'passes_all_gates': bool(gates_tuple(m))})
json.dump({'prereg': 'docs/PREREG_C3_DISULFIDE_GATE_2026-09-27.md',
           'rng_seed': 109, 'proposals': n_prop, 'gate_rejected': n_rej,
           'candidates': out,
           'scramble_control_seed111': scrambles,
           'random_panel_seed113': panel},
          open('results/c3_disulfide_screen.json', 'w'), indent=2)
n_pass = sum(1 for o in out if o['passes_all_gates'])
n_scr_pass = sum(1 for s in scrambles if s['passes_all_gates'])
n_pan_pass = sum(1 for p in panel if p['passes_all_gates'])
print(f"proposals {n_prop}, gate-rejected {n_rej}, candidates {len(out)}, passing {n_pass}")
print(f"scrambles passing {n_scr_pass}/{len(scrambles)}; random panel passing {n_pan_pass}/{len(panel)}")
for o in out[:8]:
    print(o['sequence'], round(o['score'],3), 'acp', round(o['p_acp'],3), 'sol', round(o['p_sol'],3), 'agg', round(o['p_agg'],3), 'GATES', o['passes_all_gates'])
print('elapsed', round(time.time() - t0, 1), 's')
