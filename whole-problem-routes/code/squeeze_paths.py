"""Checks for Proposition 7.23 (squeeze paths through critical points) and Observation 7.24 (one-copy compression).

Setting: f(x, y) = <y| P_a^{(x)n}(x x^T) |y> on real vectors, C = a I - (F + Omega)/4 the Choi matrix of P_a.
A joint critical point is a pair of unit vectors with P^{(x)n}(yy^T) x = f x and P^{(x)n}(xx^T) y = f y.
For a copy s, a unit w in R^3 and v in [0, 1] put G_v = I - (1 - sqrt(1 - v)) w w^T at copy s and
F(v) = f(G_v x, G_v y) / (|G_v x|^2 |G_v y|^2).  Proposition 7.23:
    F(v) - f = v^2 [ (a - 1/2) f_{n-1}(<w|x, <w|y) - f |<w|x|^2 |<w|y|^2 ] / [ (1 - v|<w|x|^2)(1 - v|<w|y|^2) ].

1. n = 3, a = 3/4 and a = 0.7: joint critical points from alternating least-eigenvector iterations (minima) and from
   minimising the squared residual of the critical equations from random starts (these are thick saddles);
   the formula is compared with direct evaluation for random w and v in {0.3, 0.8, 1}.
2. At the thick critical points found, the largest value over w of Phi(w) = f |<w|x|^2 |<w|y|^2 - (a - 1/2) f_2(...)
   (a positive value means the squeeze path in that direction lowers the normalised value).
3. Observation 7.24 (floating point): g(y) = lambda_min(P^{(x)n}(yy^T)) / |y|^2 and the best one-copy compression
   min_{s,w} g(P_{s,w} y) for n = 2 (all Schmidt coefficient vectors on a grid, using O(3) x O(3) covariance) and for
   random y at n = 3, 4.  The gap min g(Py) - g(y) should be <= 0.
Ends with PASS or FAIL (for 1 and for the sign of the gaps in 3)."""
import numpy as np
from scipy.optimize import minimize
d = 3
rng = np.random.default_rng(723)
ok = True
I9 = np.eye(9); F = np.zeros((9, 9)); Om = np.zeros((9, 9))
for i in range(3):
    for j in range(3):
        F[3*i+j, 3*j+i] = 1; Om[3*i+i, 3*j+j] = 1
def choi_power(a, m):       # C^{(x)m} with the ordering (A_1..A_m)(B_1..B_m)
    C = a*I9 - (F + Om)/4
    W = np.ones((1, 1))
    for _ in range(m): W = np.kron(W, C)
    T = W.reshape([3]*(4*m))
    perm = [2*s for s in range(m)] + [2*s+1 for s in range(m)]
    perm = perm + [2*m + p for p in perm]
    return T.transpose(perm).reshape(9**m, 9**m)
def apply_map(a, X, n):
    T = X.reshape([d] * (2 * n)).astype(float)
    for s in range(n):
        T = np.moveaxis(T, [s, n + s], [0, 1]); sh = T.shape
        Mi = T.reshape(3, 3, -1); out = np.zeros_like(Mi); tr = np.einsum('iir->r', Mi)
        for k in range(3): out[k, k] += a * tr
        out -= 0.25 * (Mi + Mi.transpose(1, 0, 2))
        T = out.reshape(sh); T = np.moveaxis(T, [0, 1], [s, n + s])
    return T.reshape(d**n, d**n)
def marg_min(v, n):
    V = v.reshape([3] * n)
    return min(np.linalg.eigvalsh(np.moveaxis(V, s, 0).reshape(3, -1) @ np.moveaxis(V, s, 0).reshape(3, -1).T)[0] for s in range(n))
# ---------------------------------------------------------------- parts 1 and 2
n = 3; D = 27; Dr = 9
for a in (0.75, 0.70):
    Cn = choi_power(a, n); Cr = choi_power(a, n - 1); T4 = Cn.reshape(D, D, D, D)
    fv = lambda x, y: np.kron(x, y) @ Cn @ np.kron(x, y)
    pts = []
    for trial in range(3):           # minima / low critical points by alternating eigenvector iteration
        y = rng.standard_normal(D); y /= np.linalg.norm(y)
        for it in range(2000):
            Wy = np.einsum('abcd,b,d->ac', T4, y, y); x = np.linalg.eigh((Wy + Wy.T)/2)[1][:, trial]
            Wx = np.einsum('abcd,a,c->bd', T4, x, x); y = np.linalg.eigh((Wx + Wx.T)/2)[1][:, trial]
        pts.append((x, y))
    def res(z):
        x = z[:D] / np.linalg.norm(z[:D]); y = z[D:] / np.linalg.norm(z[D:])
        Wy = np.einsum('abcd,b,d->ac', T4, y, y); Wx = np.einsum('abcd,a,c->bd', T4, x, x); f0 = x @ Wy @ x
        return np.sum((Wy @ x - f0 * x)**2) + np.sum((Wx @ y - f0 * y)**2)
    tries = 0
    while sum(1 for p in pts if marg_min(p[0], n) > 0.05 and marg_min(p[1], n) > 0.05) < 2 and tries < 20:
        tries += 1
        r = minimize(res, rng.standard_normal(2 * D), method='L-BFGS-B', options={'maxiter': 20000, 'gtol': 1e-14, 'ftol': 1e-20})
        if r.fun < 1e-14:
            pts.append((r.x[:D] / np.linalg.norm(r.x[:D]), r.x[D:] / np.linalg.norm(r.x[D:])))
    for (x, y) in pts:
        f0 = fv(x, y)
        Wy = np.einsum('abcd,b,d->ac', T4, y, y); Wx = np.einsum('abcd,a,c->bd', T4, x, x)
        resid = np.linalg.norm(Wy @ x - f0 * x) + np.linalg.norm(Wx @ y - f0 * y)
        X = x.reshape(3, Dr); Y = y.reshape(3, Dr)
        err = 0.0; phimax = -np.inf
        for k in range(60):
            w = rng.standard_normal(3); w /= np.linalg.norm(w)
            xi, eta = w @ X, w @ Y; al, be = xi @ xi, eta @ eta
            fp = np.kron(xi, eta) @ Cr @ np.kron(xi, eta)
            phimax = max(phimax, f0 * al * be - (a - .5) * fp)
            for v in (0.3, 0.8, 1.0):
                G = np.eye(3) - (1 - np.sqrt(1 - v)) * np.outer(w, w)
                xv, yv = (G @ X).reshape(-1), (G @ Y).reshape(-1)
                Fv = fv(xv, yv) / ((xv @ xv) * (yv @ yv))
                pred = f0 + v**2 * ((a - .5) * fp - f0 * al * be) / ((1 - v * al) * (1 - v * be))
                err = max(err, abs(Fv - pred))
        thick = min(marg_min(x, n), marg_min(y, n))
        good = err < 1e-6 * max(1.0, 1e9 * resid)
        ok &= good
        print(f"a = {a}: critical point f = {f0:+.6e} (residual {resid:.1e}, least one-copy marginal eigenvalue {thick:.3f}); "
              f"formula error {err:.1e}; max_w Phi on copy 0 (60 w) {phimax:+.3e}{'' if good else '  FAIL'}", flush=True)
# ---------------------------------------------------------------- part 3
def g(y, n, a=0.75):
    W = apply_map(a, np.outer(y, y), n); W = (W + W.T) / 2
    return np.linalg.eigvalsh(W)[0] / (y @ y)
def proj(y, s, w, n):
    w = w / np.linalg.norm(w); P = np.eye(3) - np.outer(w, w)
    V = np.moveaxis(y.reshape([d] * n), s, 0); V = np.tensordot(P, V, axes=(1, 0))
    return np.moveaxis(V, 0, s).reshape(-1)
def sph(t): return np.array([np.sin(t[0]) * np.cos(t[1]), np.sin(t[0]) * np.sin(t[1]), np.cos(t[0])])
grid = [(th, ph) for th in np.linspace(0.1, np.pi / 2, 6) for ph in np.linspace(0, 2 * np.pi, 12, endpoint=False)] + [(0.0, 0.0)]
def best_comp(y, n):
    best = np.inf
    for s in range(n):
        vals = [g(proj(y, s, sph(t), n), n) for t in grid]
        for i in np.argsort(vals)[:2]:
            r = minimize(lambda t: g(proj(y, s, sph(t), n), n), grid[i], method='Nelder-Mead',
                         options={'xatol': 1e-8, 'fatol': 1e-12, 'maxiter': 300})
            best = min(best, r.fun)
    return best
worst2 = -np.inf
for u in np.linspace(0.05, np.pi / 2 - 0.05, 12):
    for t in np.linspace(0.05, np.pi / 2 - 0.05, 12):
        s3 = np.array([np.sin(u) * np.cos(t), np.sin(u) * np.sin(t), np.cos(u)])
        y = np.diag(s3).reshape(-1)
        worst2 = max(worst2, best_comp(y, 2) - g(y, 2))
print(f"n = 2, 144 Schmidt vectors: largest gap min_(s,w) g(Py) - g(y) = {worst2:+.3e}", flush=True)
ok &= worst2 < 0
for n3, K in ((3, 8), (4, 3)):
    worst = -np.inf
    for k in range(K):
        y = rng.standard_normal(3**n3); y /= np.linalg.norm(y)
        worst = max(worst, best_comp(y, n3) - g(y, n3))
    print(f"n = {n3}, {K} random y: largest gap = {worst:+.3e}", flush=True)
    ok &= worst < 0
print("PASS" if ok else "FAIL")
