"""Adversarial search INSIDE the tangent class against the two surviving forms of the one-copy condition.
mode 'common': maximise max_{lam in [-1/2, -delta]} min_s lambda_min(Q_s - lam N_s)/||N_s||  (refutes the form Theorem 3 needs)
mode 'exists': maximise min_s max_{lam in R} lambda_min(Q_s - lam N_s)/||N_s||               (refutes SC_exists, Theorem 3')
Thickness (range and co-range local supports, every copy) is kept >= tau by a penalty.
usage: python sc_T_adv.py n mode start seed steps   (start: random | W | w3ce)"""
import numpy as np, sys, time

from common import squeeze_forms, q, tangent_range
from scipy.optimize import minimize, minimize_scalar
n = int(sys.argv[1]); mode = sys.argv[2]; start = sys.argv[3]; seed = int(sys.argv[4]); steps = int(sys.argv[5])
D = 3 ** n; rng = np.random.default_rng(seed); delta = 0.01; tau = float(sys.argv[6]) if len(sys.argv) > 6 else 0.05
P = 8 * n + 4 + 4 * D                       # U (n 2x2 complex), a (2 complex), V (2 x D complex)


def build(x):
    z = x[:P // 2] + 1j * x[P // 2:]
    U = [z[4 * t:4 * t + 4].reshape(2, 2) for t in range(n)]; a = z[4 * n:4 * n + 2]; V = z[4 * n + 2:].reshape(2, D)
    u = tangent_range(U, a, n); X = np.outer(u[0], V[0].conj()) + np.outer(u[1], V[1].conj()); return X / np.linalg.norm(X)


def thickness(X):
    out = []
    for s in range(n):
        for Y in (X @ X.conj().T, X.conj().T @ X):
            T = np.moveaxis(Y.reshape([3] * (2 * n)), [s, n + s], [0, n]).reshape(3, D // 3, 3, D // 3)
            ev = np.linalg.eigvalsh(np.einsum('aibi->ab', T)); out.append(ev[0] / ev[-1])
    return min(out)


def score(X):
    forms = [squeeze_forms(X, n, s) for s in range(n)]; sc = [np.linalg.norm(Kn, 2) for _, Kn in forms]
    if mode == "common":
        h = lambda l: min(np.linalg.eigvalsh(K - l * Kn)[0] / c for (K, Kn), c in zip(forms, sc))
        r = minimize_scalar(lambda l: -h(l), bounds=(-0.5, -delta), method='bounded', options={'xatol': 1e-9}); return -r.fun
    vals = []
    for (K, Kn), c in zip(forms, sc):
        r = minimize_scalar(lambda l: -np.linalg.eigvalsh(K - l * Kn)[0] / c, bounds=(-8, 8), method='bounded', options={'xatol': 1e-9})
        vals.append(-r.fun)
    return min(vals)


def obj(x):
    X = build(x); return -(score(X) - 1000 * max(0.0, tau - thickness(X)) ** 2)


x0 = rng.normal(size=P)
if start in ("W", "w3ce"):
    z = np.zeros(P // 2, complex)
    for t in range(n): z[4 * t:4 * t + 4] = np.eye(2).ravel()
    z[4 * n + 2:] = (rng.normal(size=2 * D) + 1j * rng.normal(size=2 * D))
    if start == "w3ce" and n == 3:
        def ket(s): v = np.zeros(D); v[9 * s[0] + 3 * s[1] + s[2]] = 1; return v
        g1 = -ket((0,0,0)) - 3*ket((0,0,1)) - 3*ket((0,1,0)) + ket((0,2,0)) + 3*ket((1,0,0)) - ket((1,0,1)) - ket((1,1,2)) - ket((1,2,1)) + 2*ket((2,1,1))
        g2 = ket((0,0,0)) - 3*ket((0,0,2)) - 3*ket((0,2,0)) - 2*ket((1,2,2)) + 3*ket((2,0,0)) + ket((2,0,2)) + ket((2,1,2)) + ket((2,2,1))
        z[4 * n + 2:] = np.concatenate([g1, g2]) + 0.05 * (rng.normal(size=2 * D) + 1j * rng.normal(size=2 * D))
    x0 = np.concatenate([z.real, z.imag]) + 0.01 * rng.normal(size=P)
X = build(x0); print(f"n={n} {mode} start {start}: score {score(X):+.4e}, thick {thickness(X):.3f}, R {q(X, n):.4f}", flush=True)
t0 = time.time(); best = [np.inf, None]


def cb(xk):
    f = obj(xk)
    if f < best[0]: best[0] = f; best[1] = xk.copy()


res = minimize(obj, x0, method='L-BFGS-B', options={'maxiter': steps, 'eps': 1e-6}, callback=cb)
X = build(best[1] if best[1] is not None else res.x)
print(f"end: score {score(X):+.4e}, thick {thickness(X):.3f}, R {q(X, n):.4f}, {time.time() - t0:.0f}s", flush=True)
