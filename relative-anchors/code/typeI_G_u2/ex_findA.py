"""Exploration (floats): locate the sup of ||Ghat_tt'||_F sigma/sqrt(s_t s_t') for u = 2 (direct elimination)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cf import TypeI

def ratio(th, ph, mu, la):
    sig = 1 - mu
    s = np.array([th * sig, ph * sig, (1 - th - ph) * sig, mu * la, mu * (1 - la)])
    w = np.sqrt(s)
    ti = TypeI(list(w), [0, 1, 2], F=float)
    M = np.array(ti.build()); U = ti.U
    tb = []
    for t in ti.T:
        tb.append(ti.idx[('R', t)])
        for h in U: tb.append(ti.idx[('Z', t, h)])
    el = [i for i in range(ti.n) if i not in tb and not (ti.labels[i][0] == 'P' and ti.labels[i][1] == U[0])]
    G = M[np.ix_(tb, tb)] - M[np.ix_(tb, el)] @ np.linalg.solve(M[np.ix_(el, el)], M[np.ix_(el, tb)])
    wU = w[3:]; S = s[3:].sum(); wh = wU / np.sqrt(S); c = np.sqrt(1 + S / sig)
    C = np.zeros((3, 3)); C[0, 0] = 1; C[1:, 1] = c * wh; C[1:, 2] = [-wh[1], wh[0]]
    Gh = C.T @ G[0:3, 3:6] @ C
    return np.linalg.norm(Gh, 'fro') * sig / np.sqrt(s[0] * s[1])

rng = np.random.default_rng(1)
best = (0, None)
for it in range(1500):
    th, ph = rng.uniform(0, 1, 2)
    if th + ph > 1: th, ph = 1 - th, 1 - ph
    th, ph = max(th, 1e-4), max(ph, 1e-4)
    mu = rng.uniform(0.001, 0.999); la = rng.choice([rng.uniform(0, 1), 1e-3, 1 - 1e-3])
    r = ratio(th, ph, mu, la)
    if r > best[0]: best = (r, (th, ph, mu, la))
print(best, flush=True)
for it in range(600):
    cand = np.clip(np.array(best[1]) + rng.normal(0, 0.03, 4), 1e-4, 1 - 1e-4)
    if cand[0] + cand[1] > 1: continue
    r = ratio(*cand)
    if r > best[0]: best = (r, tuple(cand))
print(best)
