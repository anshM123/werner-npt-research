"""Observation 7.25 (census): thick joint critical points of f_n (a = 3/4) from random starts, found by minimising the
residual of the critical equations (torch autograd through the copy-wise map, L-BFGS), and at each one the
best one-copy squeeze margin  max_{s,w} [lam - (a-1/2) f_{n-1}(xi_hat, eta_hat)]  (exact over w by Nelder-Mead).
A positive margin means that the squeeze in that direction lowers the normalised value (Proposition 7.23).
usage: python critical_census.py n seed runs   (logs: 3 21 12 and 4 31 12)"""
import sys, numpy as np, torch
from scipy.optimize import minimize
torch.set_default_dtype(torch.float64); torch.set_num_threads(4)
a = 0.75; n = int(sys.argv[1]); seed = int(sys.argv[2]); runs = int(sys.argv[3])
rng = np.random.default_rng(seed); D = 3**n
def pmap(X, m):
    # apply P_a^{(x)m} to a (3^m x 3^m) matrix X (torch)
    T = X.reshape([3] * (2 * m))
    for s in range(m):
        T = torch.movedim(T, (s, m + s), (0, 1)); sh = T.shape
        Mi = T.reshape(3, 3, -1); tr = torch.einsum('iir->r', Mi)
        out = -0.25 * (Mi + Mi.transpose(0, 1)) + a * torch.eye(3).unsqueeze(-1) * tr
        T = torch.movedim(out.reshape(sh), (0, 1), (s, m + s))
    return T.reshape(3**m, 3**m)
def resid(xr, yr):
    x = xr / xr.norm(); y = yr / yr.norm()
    Wy = pmap(torch.outer(y, y), n); Wx = pmap(torch.outer(x, x), n)
    lam = x @ Wy @ x
    return ((Wy @ x - lam * x)**2).sum() + ((Wx @ y - lam * y)**2).sum(), lam, x, y
def marg_min(v):
    V = v.reshape([3] * n)
    return min(np.linalg.eigvalsh(np.moveaxis(V, s, 0).reshape(3, -1) @ np.moveaxis(V, s, 0).reshape(3, -1).T)[0] for s in range(n))
def fprev(xi, eta):
    Wt = pmap(torch.tensor(np.outer(xi, xi)), n - 1).numpy(); return eta @ Wt @ eta
def margin(x, y, lam):
    best = -np.inf
    for s in range(n):
        Xs = np.moveaxis(x.reshape([3] * n), s, 0).reshape(3, -1); Ys = np.moveaxis(y.reshape([3] * n), s, 0).reshape(3, -1)
        def neg(t):
            w = np.array([np.sin(t[0])*np.cos(t[1]), np.sin(t[0])*np.sin(t[1]), np.cos(t[0])])
            xi, eta = w @ Xs, w @ Ys; al, be = xi @ xi, eta @ eta
            return -(lam - (a - .5) * fprev(xi, eta) / (al * be))
        for t0 in [(th, p) for th in np.linspace(0.15, 1.45, 4) for p in np.linspace(0, 6.0, 6)]:
            r = minimize(neg, t0, method='Nelder-Mead', options={'xatol': 1e-9, 'fatol': 1e-13, 'maxiter': 250})
            best = max(best, -r.fun)
    return best
out = []
for run in range(runs):
    z0 = rng.standard_normal(2 * D)
    def fr(z):
        zt = torch.tensor(z, requires_grad=True); r, lam, x, y = resid(zt[:D], zt[D:]); r.backward()
        return r.item(), zt.grad.numpy()
    r = minimize(fr, z0, jac=True, method='L-BFGS-B', options={'maxiter': 60000, 'gtol': 1e-16, 'ftol': 1e-24})
    res, lam, x, y = resid(torch.tensor(r.x[:D]), torch.tensor(r.x[D:]))
    x, y, lam = x.numpy(), y.numpy(), lam.item(); thick = min(marg_min(x), marg_min(y))
    if res.item() > 1e-14 or thick < 1e-3:
        print(f"run {run}: residual {res.item():.1e}, thickness {thick:.1e} -> skipped (value {lam:+.5f})", flush=True); continue
    m = margin(x, y, lam); out.append((lam, m, thick))
    pass
    print(f"run {run}: thick critical point lam = {lam:+.6f}, thickness {thick:.3f}, best one-copy squeeze margin {m:+.4e}"
          f"{'  <-- UPHILL EVERYWHERE' if m < 0 else ''}", flush=True)
if out:
    print(f"n = {n}: {len(out)} thick critical points, values {min(o[0] for o in out):.4f}-{max(o[0] for o in out):.4f}, "
          f"smallest margin {min(o[1] for o in out):+.4f}; product value {(a-.5)**n:.5f}")
