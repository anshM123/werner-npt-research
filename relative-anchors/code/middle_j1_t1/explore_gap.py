"""Floating-point exploration (no claim rests on it).  For the P-first reduction Phi (tlib.py) of the |T| = 1 block of
(R_inf)_1: the smallest eigenvalue of the congruence-normalised matrix Phi_n = Dg^{-1/2} Phi Dg^{-1/2}, Dg = diag(Phi),
over the U-weights at fixed sigma = x_t and fixed k, compared with its value at equal U-weights.  Also a local
minimisation (random search) of lambda_min over the U-weight simplex.  High precision (mpmath, 40 digits) avoids the
cancellations of the closed form near degenerate weights."""
import sys
import numpy as np
import mpmath as mp
from tlib import phi

mp.mp.dps = 40
rng = np.random.default_rng(3)


def lmin(xU, sig):
    """xU: U-weights (sum 1 - sig), sig = x_t; returns the smallest eigenvalue of the diagonally scaled Phi."""
    x = [mp.mpf(sig)] + [mp.mpf(v) for v in xU]
    w = [mp.sqrt(v) for v in x]
    P = phi(w, 0, one=mp.mpf(1))
    n = len(P)
    d = [mp.sqrt(P[i][i]) for i in range(n)]
    A = mp.matrix([[P[i][j] / d[i] / d[j] for j in range(n)] for i in range(n)])
    ev = mp.eigsy(A, eigvals_only=True)
    return float(min(ev))


for k in (4, 5, 7):
    u = k - 1
    for sig in (0.01, 0.3, 0.9):
        eq = lmin([(1 - sig) / u] * u, sig)
        best = (eq, 'equal')
        for trial in range(60):
            mode = trial % 3
            if mode == 0:
                y = rng.dirichlet(np.ones(u))
            elif mode == 1:
                y = rng.dirichlet(0.2 * np.ones(u))
            else:
                y = np.ones(u); y[rng.integers(u)] = rng.uniform(1e-4, 1e-2); y /= y.sum()
            y = np.maximum(y, 1e-9); y /= y.sum()
            v = lmin(list((1 - sig) * y), sig)
            if v < best[0]:
                best = (v, np.round(y, 4))
        print(f'k = {k}, sigma = {sig}: equal U-weights {eq:.4f}; smallest found {best[0]:.4f} at U-shape {best[1]}',
              flush=True)
