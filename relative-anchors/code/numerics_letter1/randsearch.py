"""Floating-point search for a negative eigenvalue of the letter-1 block of (R_inf)_j at random weights (exploration;
no claim rests on it).  Reports, per (k, j), the smallest eigenvalue above the known kernel (in units 2^Z Y) over the
samples, and the size of the eigenvalues counted as kernel."""
import sys, time
import numpy as np
from l1float import L1Block
k, nsamp, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
rng = np.random.default_rng(seed)
T0 = time.time()
for j in range(1, k - 1):
    worst = (np.inf, None); kermax = 0.0
    for it in range(nsamp):
        mode = it % 4
        if mode == 0:
            s = np.exp(rng.normal(0, 1.5, size=k))
        elif mode == 1:
            s = np.exp(rng.uniform(-6, 6, size=k))
        elif mode == 2:
            s = np.ones(k); s[rng.integers(k)] = np.exp(rng.uniform(-8, 8))
        else:
            s = np.exp(rng.normal(0, 0.3, size=k)); s[:2] *= np.exp(rng.uniform(-5, -2))
        s = s / s.max()
        B = L1Block(list(s), j)
        Y = B.build()
        ev = np.linalg.eigvalsh(Y)
        nk = B.dim(j) + (1 if j == k - 2 else 0) + (k - 1 if j == k - 1 else 0)
        kermax = max(kermax, np.abs(ev[:nk]).max())
        if ev[nk] < worst[0]:
            worst = (ev[nk], np.round(s, 4))
    print(f'k={k} j={j}: min eigenvalue above kernel {worst[0]:.3e} (kernel eigenvalues up to {kermax:.1e}); '
          f'at s = {worst[1]} [{time.time()-T0:.0f}s]', flush=True)
