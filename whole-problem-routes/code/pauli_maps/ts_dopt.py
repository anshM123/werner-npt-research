"""For lambda = (t, 0, s) at the FM boundary: is there a single product D = d^n with P = (F + D - Gamma D)/2 >= 0?
Optimise d (4 numbers) to maximise min_J P(J) (relative to the scale of F)."""
import numpy as np, itertools, sys
from scipy.optimize import minimize
from pauli_general import bell, gam, fm_ok
n = int(sys.argv[1]); rng = np.random.default_rng(int(sys.argv[2])); trials = int(sys.argv[3])
strings = [list(J) for J in itertools.product(range(4), repeat=n)]
for _ in range(trials):
    t, s = rng.uniform(0.2, 1, 2)
    d0 = np.array([t, 0.0, s])
    lo, hi = 0.0, 1.0 / max(t, s)
    for _ in range(50):
        mid = (lo + hi) / 2
        if fm_ok(mid * d0, n): lo = mid
        else: hi = mid
    l = 0.9999 * lo * d0
    if 1 - l[0] - l[2] >= 0: continue
    f = bell(l)
    def minP(d):
        gd = gam @ d
        return min(0.5 * (np.prod(f[J]) + np.prod(d[J]) - np.prod(gd[J])) for J in strings)
    best = None
    for st in [np.abs(f), np.abs(f) + 0.1, np.abs(f) * 1.2]:
        res = minimize(lambda d: -minP(d), st, method='Nelder-Mead', options={'maxiter': 3000, 'xatol': 1e-9, 'fatol': 1e-12})
        if best is None or res.fun < best.fun: best = res
    d = best.x; gd = gam @ d
    print(f"t={l[0]:.4f} s={l[2]:.4f}: f={np.round(f,4)} best min P={-best.fun:+.2e} d={np.round(d,4)} gamma d={np.round(gd,4)}")
