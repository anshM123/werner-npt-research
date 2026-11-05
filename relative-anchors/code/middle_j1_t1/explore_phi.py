"""Floating-point exploration (no claim rests on it): the P-first reduction Phi of the |T| = 1 block of (R_inf)_1.
Prints the smallest eigenvalue of Phi and of its diagonally scaled version Dg^{-1/2} Phi Dg^{-1/2} (Dg = diag(Phi))
under degenerations of the weights.  The closed form has cancellations near degenerate weights (factors such as
1 - x_c - x_u close to 0), so 40-digit arithmetic (mpmath) is used; in double precision a spurious negative eigenvalue
(-1e-4 at k = 3, w ~ (1, 0.0123, 0.0001)) appeared, and exact arithmetic shows Phi > 0 there."""
import numpy as np
import mpmath as mp
from tlib import phi

mp.mp.dps = 40
rng = np.random.default_rng(1)


def lam(w, t):
    w = [mp.mpf(float(v)) for v in w]
    nrm = mp.sqrt(sum(v * v for v in w))
    w = [v / nrm for v in w]
    P = phi(w, t, one=mp.mpf(1))
    n = len(P)
    d = [mp.sqrt(P[i][i]) for i in range(n)]
    e1 = mp.eigsy(mp.matrix(P), eigvals_only=True)
    e2 = mp.eigsy(mp.matrix([[P[i][j] / d[i] / d[j] for j in range(n)] for i in range(n)]), eigvals_only=True)
    return float(min(e1)), float(min(e2))


for k in (3, 4, 5, 6, 8):
    worst = (9, 9)
    for trial in range(60):
        mode = trial % 6
        if mode == 0:
            w = rng.uniform(0.01, 1, k)
        elif mode == 1:
            w = 3.0 ** rng.integers(-5, 5, k)
        elif mode == 2:
            w = rng.uniform(0.5, 1, k); w[rng.integers(k)] = 1e-3
        elif mode == 3:
            w = rng.uniform(0.5, 1, k); w[0] = 1e-3; w[-1] = 1e-4
        elif mode == 4:
            w = rng.uniform(0.1, 1, k); w[rng.integers(k)] = 100
        else:
            w = np.ones(k); w[rng.integers(k)] = rng.uniform(0.001, 1000)
        t = int(rng.integers(k))
        l1, l2 = lam(w, t)
        if l2 < worst[1]:
            worst = (l1, l2, np.round(np.array(w) / np.linalg.norm(w), 4), t)
    print(k, 'smallest eigenvalue of Phi and of the scaled Phi (60 samples):', worst, flush=True)

# degenerations at equal U-weights
for k in (3, 5, 8, 12):
    for sig in (1e-4, 1e-2, 0.3, 0.9, 0.999):
        w = [np.sqrt((1 - sig) / (k - 1))] * k
        w[0] = np.sqrt(sig)
        l1, l2 = lam(w, 0)
        print(f'equal U k={k} sigma={sig}: lmin Phi {l1:.4e}, scaled {l2:.4e}', flush=True)
