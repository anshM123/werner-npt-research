"""Control above the endpoint for the letter-resolved Hermitian relaxation (herm4) at n = 4: with weight w = 0.55 the
form q_n^(w)(C) = sum_J (-w)^|J| ||Tr_J C||^2 is negative on explicit Hermitian rank-2 operators, so a valid relaxation
must return a value at most that negative number.
   t = +1: C = (P_2 on one copy) (x) |x><x|^(x)3 (product):           q = (2 - 4w)(1 - w)^3
   t = -1: C = (|ab><ab| - |a'b'><a'b'|) (x) |x><x|^(x)2, a perp a', b perp b':  q = (2 - 4w)(1 - w)^2"""
import sys, os, time, numpy as np
from alg4 import build_algebras
from herm4 import Herm4, ALL_CONS
n, w = 4, float(sys.argv[1]) if len(sys.argv) > 1 else 0.55
algs = build_algebras(); t0 = time.time()
rel = Herm4(n, algs, ALL_CONS, verbose=False, w=w)
print(f"build {time.time() - t0:.0f}s (w = {w})", flush=True)
explicit = {+1.0: (2 - 4 * w) * (1 - w) ** 3, -1.0: (2 - 4 * w) * (1 - w) ** 2}
for t in (+1.0, -1.0):
    sol = rel.solve(1.0, t, max_iter=400)
    print(f"n={n} w={w} t={t:+.0f}: status {sol.status}, relaxation {sol.obj_val:+.6e}, explicit operator {explicit[t]:+.6e}  "
          f"({'OK: relaxation <= explicit' if sol.obj_val <= explicit[t] + 1e-7 else 'VIOLATION'})", flush=True)
