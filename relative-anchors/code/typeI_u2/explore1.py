"""Exploration: (C_u) as the principal block on {P_{r,h}: r, h in U} u {Q}.  Floats.
(1) the principal block is PD for random weights (u = 2..5, m = 1..3);
(2) it is PD iff the brief's (C_u) Schur complement is PD (compare smallest eigenvalues' signs; both > 0)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cf import TypeI

rng = np.random.default_rng(1)
worst = np.inf
for u in (2, 3, 4, 5):
    for m in (1, 2, 3):
        k = u + m
        T = list(range(m))
        for trial in range(60):
            spread = rng.choice([0.5, 2, 4])
            w = np.exp(rng.uniform(-spread, spread, k))
            if trial % 3 == 1:
                w[:m] *= np.exp(rng.uniform(-4, 2))
            ti = TypeI(w, T, F=float)
            M = np.array(ti.build())
            U = ti.U
            keep = [ti.idx[('P', r, h)] for r in U for h in U] + [ti.idx[('Q',)]]
            Mp = M[np.ix_(keep, keep)]
            lam = np.linalg.eigvalsh(Mp)[0]
            s = w ** 2
            rel = lam / (s[T].sum() / s.sum())
            worst = min(worst, rel)
        print(u, m, 'min lam(principal block)/(s_T/|s|) so far', worst, flush=True)
