"""Descents of a sector function mu_c with the fast sector evaluator (for n = 6).  Reports where each run ends.
usage: python search_sector_fast.py n c starts seed iters"""
import numpy as np, sys, warnings
from scipy.optimize import minimize
from sector_fast import Sector
from povm_form import povm_from
warnings.filterwarnings("ignore")
n, c, starts, seed, iters = (int(a) for a in sys.argv[1:6]); rng = np.random.default_rng(seed); S = Sector(n, c)
f = lambda x: S.mu(povm_from(x, n))
for k in range(starts):
    x0 = rng.normal(size=8 * n)
    res = minimize(f, x0, method='L-BFGS-B', options={'maxiter': iters, 'gtol': 1e-10, 'ftol': 1e-13, 'eps': 1e-7})
    A = povm_from(res.x, n); r1 = sum(1 for a in A if np.linalg.eigvalsh(a)[0] < 1e-6)
    print(f"sector descent n={n} c={c} start {k}: {f(x0):.5f} -> {res.fun:.6f}; smallest eigenvalue among A_s "
          f"{min(np.linalg.eigvalsh(a)[0] for a in A):.1e}; rank-<=1 elements {r1}/{n}; iterations {res.nit}", flush=True)
