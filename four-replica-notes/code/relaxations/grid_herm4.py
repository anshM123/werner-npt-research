"""Build the swap-reduced two-letter relaxation once at n = 4 and solve on the grid C = |u><u| + t |v><v|.
Results (value, status, residuals, min eigenvalue per block family) are printed as they come; x vectors are saved."""
import sys, time, os, numpy as np
from alg4 import build_algebras
from herm4 import Herm4, ALL_CONS
n = int(sys.argv[1]); ts = [float(a) for a in sys.argv[2].split(',')]
cons = sys.argv[3].split(',') if len(sys.argv) > 3 else ALL_CONS
tag = sys.argv[4] if len(sys.argv) > 4 else 'all'
algs = build_algebras()
t0 = time.time()
rel = Herm4(n, algs, cons)
print(f"build done in {time.time() - t0:.0f}s; parameters {rel.N}, identity rank {rel.Q.shape[0]}, svec rows {rel.nrows}, "
      f"stored entries {rel.nnz / 1e6:.1f}M", flush=True)
for t in ts:
    sol = rel.solve(1.0, t, max_iter=400)
    x = np.array(sol.x); mins, mom, eqr, obj = rel.evaluate(x, 1.0, t)
    print(f"n={n} t={t:+.2f}: status {sol.status}, iters {sol.iterations}, q relaxation {sol.obj_val:.8e} (c.x {obj:.8e}), "
          f"moment res {mom:.1e}, identity res {eqr:.1e}, min eig {', '.join(f'{k}:{v:.1e}' for k, v in mins.items())}; "
          f"solve {rel.solve_time:.0f}s", flush=True)
    np.savez(f"herm4_n{n}_{tag}_t{t:+.2f}.npz", x=x, z=np.array(sol.z), s=np.array(sol.s), obj=sol.obj_val)
