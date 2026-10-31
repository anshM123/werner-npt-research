"""Section 3 (numerical observations, local optimisation only; not certified).
   (a) min of q_n(1/2, C)/||C||^2 over rank-2 C in End(Sym^n C^3), n = 2..NMAX;
   (b) along the Werner family t in {1/3, 0.4, 0.45, 0.48, 0.49, 0.5}: the min over ALL rank-2 bosonic C versus the
       min over qubit-supported ones (range and co-range in Sym^n(span{e0, e1})), n = 2..4."""
import sys
import numpy as np
from bosonic import basis, min_rank2

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6
rng = np.random.default_rng(2)
print("(a) min over rank-2 bosonic C of q_n/||C||^2 at t = 1/2")
for n in range(2, NMAX + 1):
    print(f"  n={n:2d} (dim Sym^n = {len(basis(n)):3d}): {min_rank2(n, 0.5, 12, rng=rng):+.2e}", flush=True)
print("(b) Werner family: all rank-2 bosonic C vs qubit-supported C")
for n in (2, 3, 4):
    qub = np.array([b[2] == 0 for b in basis(n)])
    for t in (1 / 3, 0.4, 0.45, 0.48, 0.49, 0.5):
        a = min_rank2(n, t, 12, rng=rng); b = min_rank2(n, t, 6, mask=qub, rng=rng)
        print(f"  n={n} t={t:.4f}: all {a:+.6f}   qubit-supported {b:+.6f}   n 2^(1-n) (1-2t) = "
              f"{n * 2 ** (1 - n) * (1 - 2 * t):.6f}", flush=True)
