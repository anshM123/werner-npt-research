"""Observation 7.25, second search: thick joint critical points of f_3(x, y) = <y|P_{3/4}^{(x)3}(xx^T)|y> and their
best one-copy squeeze margin m = max over copies s and unit w of [ f - (a - 1/2) f_2(xi_hat, eta_hat) ], where
xi = <w|_s x, eta = <w|_s y and hats mean normalised.  By Proposition 7.23, m > 0 means the squeeze in that direction
lowers the normalised value along its whole path; m < 0 at every s, w would mean uphill along every one-copy squeeze.
Stage 1 pulls the value toward a target with a thickness floor; stage 2 solves the critical equations alone; a point is
kept if the residual is below 1e-12 and both vectors have least one-copy marginal eigenvalue above 1e-3.
usage: python thick_critical_points.py seed targets reps   (the log used: 5 0.05,0.1,0.15 3)"""
import sys, numpy as np
from scipy.optimize import minimize
a = 0.75; n = 3; D = 27; Dr = 9
rng = np.random.default_rng(int(sys.argv[1])); targets = [float(t) for t in sys.argv[2].split(',')]; reps = int(sys.argv[3])
I9 = np.eye(9); F = np.zeros((9, 9)); Om = np.zeros((9, 9))
for i in range(3):
    for j in range(3):
        F[3*i+j, 3*j+i] = 1; Om[3*i+i, 3*j+j] = 1
C = a*I9 - (F + Om)/4
def choi_power(m):
    W = np.ones((1, 1))
    for _ in range(m): W = np.kron(W, C)
    T = W.reshape([3]*(4*m))
    perm = [2*s for s in range(m)] + [2*s+1 for s in range(m)]
    perm = perm + [2*m + p for p in perm]
    return T.transpose(perm).reshape(9**m, 9**m)
T4 = choi_power(n).reshape(D, D, D, D); T2 = choi_power(n - 1).reshape(Dr, Dr, Dr, Dr)
def unpack(z): return z[:D] / np.linalg.norm(z[:D]), z[D:] / np.linalg.norm(z[D:])
def marg_eigs(v):
    V = v.reshape(3, 3, 3); return [np.linalg.eigvalsh(np.moveaxis(V, s, 0).reshape(3, -1) @ np.moveaxis(V, s, 0).reshape(3, -1).T)[0] for s in range(3)]
def resid(z):
    x, y = unpack(z)
    Wy = np.einsum('abcd,b,d->ac', T4, y, y); Wx = np.einsum('abcd,a,c->bd', T4, x, x); f0 = x @ Wy @ x
    return np.sum((Wy @ x - f0 * x)**2) + np.sum((Wx @ y - f0 * y)**2), f0
def stage1(z, target, floor=0.03):
    r, f0 = resid(z); x, y = unpack(z)
    pen = sum(max(0, floor - e)**2 for e in marg_eigs(x) + marg_eigs(y))
    return r + 0.3 * (f0 - target)**2 + 50 * pen
def f2(xi, eta): return np.einsum('abcd,a,b,c,d->', T2, xi, eta, xi, eta)
def contraction(v, s, w):
    V = np.moveaxis(v.reshape(3, 3, 3), s, 0); return np.tensordot(w, V, axes=(0, 0)).reshape(-1)
def sph(t): return np.array([np.sin(t[0]) * np.cos(t[1]), np.sin(t[0]) * np.sin(t[1]), np.cos(t[0])])
def best_margin(x, y, f0):
    best = -np.inf
    for s in range(3):
        def neg(t):
            w = sph(t); xi = contraction(x, s, w); eta = contraction(y, s, w); al, be = xi @ xi, eta @ eta
            if al < 1e-12 or be < 1e-12: return 0.0
            return -(f0 - (a - .5) * f2(xi, eta) / (al * be))
        grid = [(th, ph) for th in np.linspace(0.1, np.pi / 2, 7) for ph in np.linspace(0, 2 * np.pi, 14, endpoint=False)] + [(0.0, 0.0)]
        vals = [neg(t) for t in grid]
        for i in np.argsort(vals)[:3]:
            r = minimize(neg, grid[i], method='Nelder-Mead', options={'xatol': 1e-10, 'fatol': 1e-14, 'maxiter': 400})
            best = max(best, -r.fun)
    return best
kept = []
for target in targets:
    for rep in range(reps):
        z = rng.standard_normal(2 * D)
        r1 = minimize(stage1, z, args=(target,), method='L-BFGS-B', options={'maxiter': 6000})
        r2 = minimize(lambda q: resid(q)[0], r1.x, method='L-BFGS-B', options={'maxiter': 20000, 'gtol': 1e-15, 'ftol': 1e-22})
        res, f0 = resid(r2.x); x, y = unpack(r2.x); thick = min(marg_eigs(x) + marg_eigs(y))
        if res > 1e-12 or thick < 1e-3:
            print(f"target {target} run {rep}: discarded (residual {res:.1e}, least marginal eigenvalue {thick:.1e}, f {f0:+.4f})", flush=True)
            continue
        m = best_margin(x, y, f0); kept.append((f0, m))
        print(f"target {target} run {rep}: thick critical point f = {f0:+.6f} (residual {res:.1e}, least marginal eigenvalue "
              f"{thick:.3f}); best one-copy squeeze margin {m:+.3e}", flush=True)
if kept:
    print(f"{len(kept)} thick critical points, values {min(k[0] for k in kept):.3f} to {max(k[0] for k in kept):.3f}; margins "
          f"{min(k[1] for k in kept):+.3f} to {max(k[1] for k in kept):+.3f}; uphill along every one-copy squeeze at "
          f"{sum(1 for k in kept if k[1] <= 0)} of them")
