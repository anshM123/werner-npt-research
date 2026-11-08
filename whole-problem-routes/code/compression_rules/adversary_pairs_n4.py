"""Adversary for the two-copy thinning statement at n = 4: maximise
   Gap(y) = min over pairs (s,t) and plane normals (w_s, w_t) of g(G_s G_t y) - g(y)
by Danskin ascent in y (inner minimiser warm-started and re-optimised every step).
g(y) = lambda_min(P^{(x)4}(yy^T))/|y|^2;  grad g(y) = 2 (W_z y - g y)/|y|^2 with z the bottom eigenvector of W_y."""
import sys, time, math, itertools
import numpy as np
from scipy.optimize import minimize
n = 4; D = 81
def Pn(X):
    T = X.reshape((3,) * (2 * n))
    for s in range(n):
        Tm = np.moveaxis(T, [s, n + s], [0, 1])
        tr = np.einsum('ii...->...', Tm)
        out = -0.25 * (Tm + np.swapaxes(Tm, 0, 1))
        for i in range(3):
            out[i, i] += 0.75 * tr
        T = np.moveaxis(out, [0, 1], [s, n + s])
    return T.reshape(D, D)
def g_and_grad(y):
    W = Pn(np.outer(y, y)); ev, V = np.linalg.eigh(W); z = V[:, 0]; nrm = y @ y
    gval = ev[0] / nrm
    grad = 2 * (Pn(np.outer(z, z)) @ y - gval * y) / nrm
    return gval, grad
def unit(th):
    return np.array([np.sin(th[0]) * np.cos(th[1]), np.sin(th[0]) * np.sin(th[1]), np.cos(th[0])])
def proj(s, t, ws, wt):
    Ps = np.eye(3) - np.outer(ws, ws); Pt = np.eye(3) - np.outer(wt, wt)
    ops = [np.eye(3)] * n; ops = list(ops); ops[s] = Ps; ops[t] = Pt
    M = ops[0]
    for O in ops[1:]:
        M = np.kron(M, O)
    return M
PAIRS = list(itertools.combinations(range(n), 2))
def inner(y, warm, rng):
    best = (np.inf, None)
    for (s, t) in PAIRS:
        starts = [warm.get((s, t))] if warm.get((s, t)) is not None else []
        starts += [rng.uniform(0, np.pi, 4) * [1, 2, 1, 2] for _ in range(2)]
        for th0 in starts:
            f = lambda th: (lambda z: g_and_grad(z)[0] if z @ z > 1e-12 else 1.0)(proj(s, t, unit(th[:2]), unit(th[2:])) @ y)
            r = minimize(f, th0, method='Nelder-Mead', options={'xatol': 1e-7, 'fatol': 1e-12, 'maxfev': 500})
            if r.fun < best[0]:
                best = (r.fun, (s, t, r.x))
            if warm.get((s, t)) is None or r.fun < f(warm[(s, t)]):
                warm[(s, t)] = r.x
    return best
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
def omega_pairs(pairs):
    v = np.zeros([3] * n)
    for idx in itertools.product(range(3), repeat=n):
        if all(idx[s] == idx[t] for s, t in pairs):
            v[idx] = 1
    return v.reshape(-1)
p1, p2, p3 = (omega_pairs(pr) for pr in ([(0, 1), (2, 3)], [(0, 2), (1, 3)], [(0, 3), (1, 2)]))
c = (math.sqrt(6805) - 71) / 36
partner = p2 - c * (p1 + p3); partner /= np.linalg.norm(partner)
starts = [('partner + noise', partner + 0.05 * rng.standard_normal(D) / 9), ('random', rng.standard_normal(D)),
          ('random', rng.standard_normal(D)), ('partner + noise', partner + 0.15 * rng.standard_normal(D) / 9)]
t0 = time.time(); overall = -np.inf
for name, y in starts:
    y = y / np.linalg.norm(y); warm = {}; step = 0.05
    for it in range(60):
        gy, grad_y = g_and_grad(y)
        gc, (s, t, th) = inner(y, warm, rng)
        gap = gc - gy
        P = proj(s, t, unit(th[:2]), unit(th[2:]))
        gz, grad_z = g_and_grad(P @ y)
        G = P.T @ grad_z - grad_y           # ascent direction for the gap with the inner minimiser fixed
        G -= (G @ y) * y
        y_new = y + step * G / (np.linalg.norm(G) + 1e-15); y_new /= np.linalg.norm(y_new)
        y = y_new
        if it % 10 == 0 or it == 59:
            print(f'  ({name}) it {it}: g(y) = {gy * 256:.4f}, best two-copy compression gap {gap * 256:+.5f} x 4^-4 at pair ({s},{t}) [{time.time() - t0:.0f}s]', flush=True)
    overall = max(overall, gap)
    pass  # final vectors not kept
print(f'largest final gap {overall * 256:+.5f} x 4^-4 (positive would refute the two-copy statement)')
