"""Floating-point test (no claim rests on it): m = 1, u = 1..8, v = (1, beta, ..., beta) with
beta = 3 (2 s_t + S) / (4 s_t (s_t + S)) (the vector of check_mvec.py).  Checks, over random weights with spreads up to
e^{+-7} and shrunk families: G_z is a Z-matrix and G_z v > 0 entrywise.  Prints the worst normalised margins."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from numG import Gmat
rng = np.random.default_rng(23)
allok = True
for u in range(1, 9):
    worst_off = -np.inf
    worst_row = np.inf
    N = 150 if u <= 5 else 60
    for it in range(N):
        st = np.exp(rng.uniform(-7, 7))
        sU = np.exp(rng.uniform(-7, 7, u))
        if it % 3 == 1:
            sU = sU * np.exp(-rng.uniform(0, 9))
        if it % 3 == 2:
            sU[0] *= np.exp(rng.uniform(0, 9))
        G = Gmat([st], list(sU))
        w = np.sqrt(np.array([st] + list(sU)))
        Dg = np.diag([1.0] + [w[0] * w[1 + h] for h in range(u)])
        Gz = Dg @ G @ Dg
        S = sU.sum()
        beta = 3 * (2 * st + S) / (4 * st * (st + S))
        v = np.array([1.0] + [beta] * u)
        sc = np.sqrt(np.outer(np.diag(Gz), np.diag(Gz)))
        off = (Gz / sc)[~np.eye(u + 1, dtype=bool)]
        worst_off = max(worst_off, off.max())
        rows = (Gz @ v) / (np.diag(Gz) * v)        # (G_z v)_i / (G_ii v_i): in (0, 1] when the row is dominated
        worst_row = min(worst_row, rows.min())
    ok = worst_off < 0 and worst_row > 0
    allok &= ok
    print(f'u = {u} ({N} samples): max normalised off-diagonal entry {worst_off:.2e}; min (G_z v)_i/(G_ii v_i) = {worst_row:.3e}  {"ok" if ok else "FAILS"}', flush=True)
print('PASS' if allok else 'FAIL')
