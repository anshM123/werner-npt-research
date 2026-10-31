"""Section 3 (numerical observations, local optimisation only; not certified).
   (a) min of q_n(1/2, C)/||C||^2 over rank-2 C in End(Sym^n C^3) for n = 2..NA;
   (b) Werner family: min over all rank-2 bosonic C vs qubit-supported C, for n = 2..NB."""
import sys
import numpy as np
from bosonic import basis, min_rank2

NA = int(sys.argv[1]) if len(sys.argv) > 1 else 7
NB = int(sys.argv[2]) if len(sys.argv) > 2 else 5
rng = np.random.default_rng(2)
print("(a) min over rank-2 bosonic C of q_n/||C||^2 at t = 1/2")
for n in range(2, NA + 1):
    print(f"  n={n:2d} (dim Sym^n = {len(basis(n)):3d}): {min_rank2(n, 0.5, 12, rng=rng):+.2e}", flush=True)
print("(b) Werner family: all rank-2 bosonic C vs qubit-supported C; (1-t)^n - t^n is the value at the GHZ operator")
for n in range(2, NB + 1):
    qub = np.array([b[2] == 0 for b in basis(n)])
    for t in (1 / 3, 0.4, 0.45, 0.48, 0.49, 0.5):
        a = min_rank2(n, t, 12, rng=rng); b = min_rank2(n, t, 6, mask=qub, rng=rng)
        print(f"  n={n} t={t:.4f}: all {a:+.6f}   qubit-supported {b:+.6f}   (1-t)^n - t^n = {(1 - t) ** n - t ** n:.6f}",
              flush=True)
