"""Adversarial evidence for Conjecture 12 (README): minimise lambda_min(Y_c(sum_t w_t |1_t>)) over the weights with
Nelder-Mead from random, one-dominant and two-dominant starts.  Usage: python adversarial_relative.py n [seed].
The kernel conj(p) (x) H_c makes 0 the best possible value; a clearly negative value would refute the conjecture."""
import numpy as np, sys
from scipy.optimize import minimize
from relative_anchor import Relative
n = int(sys.argv[1]); rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
worst_all = 9
for c in range(1, n + 1):
    R = Relative(n, c)
    f = lambda x: np.linalg.eigvalsh(R.Y(R.pvec(np.exp(x)))[0])[0]
    best = 9
    for k in range(5):
        x0 = rng.normal(size=n) * (2 if k % 2 else 0.5)
        if k == 3: x0 = np.array([0.0] + [-4.0] * (n - 1)) + 0.1 * rng.normal(size=n)
        if k == 4: x0 = np.array([0.0, 0.0] + [-3.0] * (n - 2)) + 0.1 * rng.normal(size=n)
        r = minimize(f, x0, method='Nelder-Mead', options=dict(maxiter=200 * n, xatol=1e-6, fatol=1e-14))
        best = min(best, r.fun)
    worst_all = min(worst_all, best)
    print(f"n={n} c={c}: adversarial min lambda_min(Y_c) = {best:+.3e}", flush=True)
print("PASS" if worst_all > -1e-12 else "FAIL")
