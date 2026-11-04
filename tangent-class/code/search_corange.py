"""W-type range fixed (U_t = I, a = 0); maximise nu_0 over the co-range only.  Counts starts that reach nu_0 > 0 with
thickness (range and co-range, every copy) >= 0.05.  usage: python corange_search.py n starts seed"""
import numpy as np, sys

from common import squeeze_forms, tangent_range
from scipy.optimize import minimize, minimize_scalar
n, starts, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); D = 3 ** n; rng = np.random.default_rng(seed)
u = tangent_range([np.eye(2)] * n, np.zeros(2), n)
def build(x):
    v = (x[:2 * D] + 1j * x[2 * D:]).reshape(2, D); X = np.outer(u[0], v[0].conj()) + np.outer(u[1], v[1].conj()); return X / np.linalg.norm(X)
def thick(X):
    out = []
    for s in range(n):
        for Y in (X @ X.conj().T, X.conj().T @ X):
            T = np.moveaxis(Y.reshape([3] * (2 * n)), [s, n + s], [0, n]).reshape(3, D // 3, 3, D // 3)
            ev = np.linalg.eigvalsh(np.einsum('aibi->ab', T)); out.append(ev[0] / ev[-1])
    return min(out)
def nu0(X):
    K, Kn = squeeze_forms(X, n, 0); c = np.linalg.norm(Kn, 2)
    r = minimize_scalar(lambda l: -np.linalg.eigvalsh(K - l * Kn)[0] / c, bounds=(-8, 8), method='bounded', options={'xatol': 1e-10}); return -r.fun
hits = 0
for k in range(starts):
    x0 = rng.normal(size=4 * D)
    res = minimize(lambda x: -(nu0(build(x)) - 1000 * max(0, 0.05 - thick(build(x))) ** 2), x0, method='L-BFGS-B', options={'maxiter': 60, 'eps': 1e-6})
    X = build(res.x); v, th = nu0(X), thick(X); hits += (v > 0 and th >= 0.05 - 1e-9)
    print(f"n={n} start {k}: nu_0 = {v:+.4e}, thickness {th:.3f}", flush=True)
print(f"n={n}: {hits}/{starts} starts reach nu_0 > 0 at thickness >= 0.05")
