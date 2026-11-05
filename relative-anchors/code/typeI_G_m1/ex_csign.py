"""Exploration (floating point): sign pattern of the eliminated block C (gauge x_t = 0, scaled coordinates
x_rh = X_rh/(w_r w_h)) and of the couplings, over random weights."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cf import TypeI
np.set_printoptions(linewidth=220, precision=3, suppress=True)

def blocks(st, sU, gauge='t'):
    w = np.sqrt(np.array([st] + list(sU), dtype=float))
    ti = TypeI(w, [0], F=float)
    M = np.array(ti.build())
    U = ti.U
    sc = np.array([w[l[1]] * w[l[2]] if l[0] in ('P', 'Z') else 1.0 for l in ti.labels])
    Ms = M * np.outer(sc, sc)
    if gauge == 't':
        el = [ti.idx[('P', r, h)] for r in U for h in U] + [ti.idx[('Q',)]]
    else:
        el = [ti.idx[('P', r, h)] for r in [0] + U[1:] for h in U] + [ti.idx[('Q',)]]
    tb = [ti.idx[('R', 0)]] + [ti.idx[('Z', 0, h)] for h in U]
    return Ms[np.ix_(el, el)], Ms[np.ix_(el, tb)], [ti.labels[i] for i in el]

if __name__ == "__main__":
    rng = np.random.default_rng(2)
    u = 3
    pats = []
    for it in range(300):
        st = np.exp(rng.uniform(-4, 4)); sU = np.exp(rng.uniform(-4, 4, u))
        C, B, lab = blocks(st, sU)
        pats.append((np.sign(C), np.sign(B)))
    SC = np.array([p[0] for p in pats]); SB = np.array([p[1] for p in pats])
    mixedC = (SC.max(0) > 0) & (SC.min(0) < 0)
    mixedB = (SB.max(0) > 0) & (SB.min(0) < 0)
    print('labels', lab)
    print('C sign (constant entries: +1/-1, mixed: 9, zero: 0):')
    P = np.where(mixedC, 9, SC[0]); print(P.astype(int))
    print('B sign:'); print(np.where(mixedB, 9, SB[0]).astype(int))
