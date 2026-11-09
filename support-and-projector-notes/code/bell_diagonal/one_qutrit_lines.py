"""Adversarial test: minimise (rank-two minimum of K_h) over h with all 12 line conditions satisfied and sum h = 9.
Random-direction search on the feasible set plus local perturbation.  usage: python one_qutrit_lines.py iters seed"""
import sys, itertools
import numpy as np
from scipy.optimize import minimize
iters, seed = int(sys.argv[1]), int(sys.argv[2])
rng = np.random.default_rng(seed)
om = np.exp(2j * np.pi / 3)
X = np.roll(np.eye(3), 1, axis=0); Zm = np.diag([1, om, om ** 2])
pts = list(itertools.product(range(3), repeat=2))
Dm = [np.linalg.matrix_power(X, z[0]) @ np.linalg.matrix_power(Zm, z[1]) for z in pts]
lines = []
for dirn in [(0, 1), (1, 0), (1, 1), (1, 2)]:
    for base in pts:
        l = frozenset(((base[0] + c * dirn[0]) % 3, (base[1] + c * dirn[1]) % 3) for c in range(3))
        if l not in lines: lines.append(l)
L_idx = [[pts.index(z) for z in l] for l in lines]
def line_margin(h):
    return min(h[i].sum() - np.linalg.norm(h[i]) for i in L_idx)
# Pauli coefficient map: c_z = Tr(D_z^* C), vectorised
Mz = np.array([d.conj().T.reshape(-1) for d in Dm])     # 9 x 9 : c = Mz @ vec(C^T)?  use einsum below
def pauli_dist(C):
    c = np.array([np.trace(d.conj().T @ C) for d in Dm])
    return np.abs(c) ** 2 / (3 * np.vdot(C, C).real)
def rank2_min(h, restarts=24):
    best = np.inf; arg = None
    for r in range(restarts):
        x0 = rng.normal(size=24)
        def f(x):
            A = (x[:6] + 1j * x[6:12]).reshape(3, 2); B = (x[12:18] + 1j * x[18:24]).reshape(3, 2)
            return h @ pauli_dist(A @ B.conj().T)
        res = minimize(f, x0, method='BFGS')
        if res.fun < best: best, arg = res.fun, res.x
    return best, arg
def project_feasible(h):
    # scale the deviation from the uniform vector until line conditions hold
    u = np.ones(9); d = h - h.mean()
    lo, hi = 0.0, 1.0
    if line_margin(u + d) >= 0: return u + d
    for _ in range(50):
        mid = (lo + hi) / 2
        if line_margin(u + mid * d) >= 0: lo = mid
        else: hi = mid
    return u + lo * d
best = (np.inf, None)
h = project_feasible(rng.normal(size=9) * 3)
cur, _ = rank2_min(h)
for it in range(iters):
    h2 = project_feasible(h + rng.normal(size=9) * (0.3 if it % 3 else 1.5))
    m, _ = rank2_min(h2, restarts=12)
    if m < cur:
        h, cur = h2, m
    if cur < best[0]:
        best = (cur, h.copy())
print(f"adversarial search ({iters} steps): smallest rank-two minimum with all line conditions satisfied = {best[0]:.3e}")
print("  h =", np.round(best[1], 5).tolist(), " line margin", f"{line_margin(best[1]):.2e}")
m, x = rank2_min(best[1], restarts=60)
print(f"  re-check with 60 restarts: {m:.3e}")
