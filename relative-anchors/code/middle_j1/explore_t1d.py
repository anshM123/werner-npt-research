"""Exploration: |T| = 1 block, rescaled.  Pi[r,r'] < 0 for r != r'.  Sign and size of the correction
b_r^T C^{-1} b_r' (C = the (R, Q) block).  Also: the ratio S_P[r,r'] / Pi[r,r']."""
import numpy as np
from j1lib import J1, block_states, mu

rng = np.random.default_rng(8)
neg = 0; tot = 0; rmin = 1e9; rmax = -1e9; cmin = 1e9
for trial in range(400):
    k = int(rng.integers(3, 10))
    x = rng.dirichlet(np.ones(k) * rng.choice([0.2, 1, 4])) + 1e-12
    w = np.sqrt(x)
    t = int(rng.integers(0, k))
    J = J1(list(w), F=float)
    st = block_states(k, (t,))
    M = np.array(J.matrix(st), dtype=float)
    m = np.array([mu(z, list(w)) for z in st])
    Mr = m[:, None] * M * m[None, :]
    P = [i for i, (a, y) in enumerate(st) if a is not None and a[0] == 1]
    E = [i for i in range(len(st)) if i not in P]
    A = Mr[np.ix_(P, P)]; B = Mr[np.ix_(P, E)]; C = Mr[np.ix_(E, E)]
    corr = B @ np.linalg.solve(C, B.T)
    S = A - corr
    for i in range(k):
        for j in range(k):
            if i != j:
                tot += 1
                neg += corr[i, j] < 0
                rmin = min(rmin, S[i, j] / A[i, j]); rmax = max(rmax, S[i, j] / A[i, j])
                cmin = min(cmin, corr[i, j] / abs(A[i, j]))
print('pairs', tot, 'with negative correlation b_r^T C^-1 b_r\':', neg)
print('ratio S_P[r,r\']/Pi[r,r\'] in', rmin, rmax, '; min corr/|Pi|', cmin)
