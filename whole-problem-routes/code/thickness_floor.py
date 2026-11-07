"""Observation 7.27(a): how much does thickness cost?  For n = 2..5 (a = 3/4) minimise f_n(x, y) over unit x, y subject to every one-copy marginal
of x and of y having least eigenvalue >= tau (penalty, then a check), for tau in a grid; compare with 4^-n and fit
min f ~ c(tau)^n.  Values are upper bounds for the constrained minima (local searches, a few starts).
usage: python thickness_floor.py seed starts   (the log used: 1 6)"""
import sys, numpy as np, torch
from scipy.optimize import minimize
torch.set_default_dtype(torch.float64); torch.set_num_threads(6)
a = 0.75
def pmap(X, m):
    T = X.reshape([3] * (2 * m))
    for s in range(m):
        T = torch.movedim(T, (s, m + s), (0, 1)); sh = T.shape
        Mi = T.reshape(3, 3, -1); tr = torch.einsum('iir->r', Mi)
        out = -0.25 * (Mi + Mi.transpose(0, 1)) + a * torch.eye(3).unsqueeze(-1) * tr
        T = torch.movedim(out.reshape(sh), (0, 1), (s, m + s))
    return T.reshape(3**m, 3**m)
def margmin(v, n):
    out = []
    for s in range(n):
        V = torch.movedim(v.reshape([3] * n), s, 0).reshape(3, -1); out.append(torch.linalg.eigvalsh(V @ V.T)[0])
    return torch.stack(out)
rng = np.random.default_rng(int(sys.argv[1])); starts = int(sys.argv[2])
for n in (2, 3, 4, 5):
    D = 3**n
    for tau in (0.02, 0.05, 0.1, 0.2, 0.3):
        best = np.inf
        for st in range(starts):
            def obj(z):
                zt = torch.tensor(z, requires_grad=True); x = zt[:D] / zt[:D].norm(); y = zt[D:] / zt[D:].norm()
                f = y @ pmap(torch.outer(x, x), n) @ y
                pen = (torch.relu(tau - margmin(x, n))**2).sum() + (torch.relu(tau - margmin(y, n))**2).sum()
                val = f + 200 * pen; val.backward()
                return val.item(), zt.grad.numpy()
            r = minimize(obj, rng.standard_normal(2 * D), jac=True, method='L-BFGS-B', options={'maxiter': 4000})
            zt = torch.tensor(r.x); x = zt[:D] / zt[:D].norm(); y = zt[D:] / zt[D:].norm()
            f = (y @ pmap(torch.outer(x, x), n) @ y).item(); th = min(margmin(x, n).min().item(), margmin(y, n).min().item())
            if th >= tau - 2e-3: best = min(best, f)
        print(f"n = {n}, tau = {tau}: least value found {best:.6f}; per-copy factor best^(1/n) = {best**(1/n):.4f}; "
              f"product value 4^-n = {4.0**-n:.6f}", flush=True)
