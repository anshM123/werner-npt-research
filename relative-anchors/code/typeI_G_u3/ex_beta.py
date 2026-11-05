"""Exploration (floating point): m = 1, which simple beta puts v = (1, beta, ..., beta) in the window
max_h |G_z[h, rho]| / r_h < beta < G_z[rho, rho] / |sum_h G_z[rho, h]|  (r_h = z-block row sums)?
Prints, for candidate beta(s_t, S), the worst ratios beta/hi (< 1 needed) and lo/beta (< 1 needed) over samples."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from numG import Gmat
rng = np.random.default_rng(11)
cands = {
    '1/s_t': lambda st, S, sU: 1 / st,
    '1/S': lambda st, S, sU: 1 / S,
    '1/(s_t+S)': lambda st, S, sU: 1 / (st + S),
    '2/(s_t+S)': lambda st, S, sU: 2 / (st + S),
    '1/(2s_t+S)': lambda st, S, sU: 1 / (2 * st + S),
    '2/(2s_t+S)': lambda st, S, sU: 2 / (2 * st + S),
    '1/(s_t+S/2)': lambda st, S, sU: 1 / (st + S / 2),
    '4/(2s_t+S)': lambda st, S, sU: 4 / (2 * st + S),
    '1/(s_t+2S)': lambda st, S, sU: 1 / (st + 2 * S),
    '3/(2(s_t+S))': lambda st, S, sU: 1.5 / (st + S),
}
worst = {k: [0.0, 0.0] for k in cands}
for u in range(1, 8):
    for it in range(100):
        st = np.exp(rng.uniform(-7, 7))
        sU = np.exp(rng.uniform(-7, 7, u))
        if it % 3 == 1:
            sU = sU * np.exp(-rng.uniform(0, 8))
        G = Gmat([st], list(sU))
        w = np.sqrt(np.array([st] + list(sU)))
        Dg = np.diag([1.0] + [w[0] * w[1 + h] for h in range(u)])
        Gz = Dg @ G @ Dg
        r = Gz[1:, 1:].sum(axis=1)
        lo = max(abs(Gz[1 + h, 0]) / r[h] for h in range(u))
        hi = Gz[0, 0] / abs(Gz[0, 1:].sum())
        S = sU.sum()
        for k, f in cands.items():
            b = f(st, S, sU)
            worst[k][0] = max(worst[k][0], b / hi)
            worst[k][1] = max(worst[k][1], lo / b)
for k, (a, b) in worst.items():
    print(f'beta = {k:14s}: max beta/hi = {a:.3f}, max lo/beta = {b:.3f}  {"WORKS" if a < 1 and b < 1 else ""}')
