"""Checks for Proposition 7.23(ii) (trace of the one-copy squeeze Hessian at critical points) and the ratios of
Observation 7.25.  At a joint critical point (x, y) of f_n with value lam, let Q_s(H) be the second derivative at t = 0
of F(t) = f(e^{tH}x, e^{tH}y)/(|e^{tH}x|^2 |e^{tH}y|^2), H real symmetric acting at copy s.  Claimed:
    trace of Q_s over traceless H = -4 [ lam (4a - 1 + 2 Tr rho^x_s rho^y_s) - (a - 1/2)(1 + 4a) tau_s ],
tau_s = <rho^x_{-s}, P_a^{(x)(n-1)}(rho^y_{-s})>.  Part 1 compares with the Hessian computed from its closed form
Q(H) = 2f(Hx,y) + 2f(x,Hy) + 4<Hx(x)y|C|x(x)Hy> + 4<Hx(x)Hy|C|x(x)y> - 2lam(|Hx|^2+|Hy|^2) - 8lam<x|H|x><y|H|y>
and with finite differences, at critical points of several kinds (products, the plane GHZ pair, thick saddles) for
n = 3 and a = 3/4, 0.7.  Part 2: for thick critical points from random starts (n = 3, 4; a = 3/4) the ratio
r_s = (a - 1/2)(1 + 4a) tau_s / [lam (4a - 1 + 2 Tr rho rho)] at every copy (r_s < 1: negative trace, not a local minimum).
Ends with PASS or FAIL (part 1 agreement, and r_s < 1 at every copy of every thick critical point in part 2)."""
import numpy as np, torch
from scipy.optimize import minimize
torch.set_default_dtype(torch.float64); torch.set_num_threads(4)
ok = True
def make(a):
    def pmap(X, m):
        T = X.reshape([3] * (2 * m))
        for s in range(m):
            T = torch.movedim(T, (s, m + s), (0, 1)); sh = T.shape
            Mi = T.reshape(3, 3, -1); tr = torch.einsum('iir->r', Mi)
            out = -0.25 * (Mi + Mi.transpose(0, 1)) + a * torch.eye(3).unsqueeze(-1) * tr
            T = torch.movedim(out.reshape(sh), (0, 1), (s, m + s))
        return T.reshape(3**m, 3**m)
    return pmap
def critical(pmap, n, rng, runs):
    D = 3**n; pts = []
    for run in range(runs):
        def fr(z):
            zt = torch.tensor(z, requires_grad=True); x = zt[:D] / zt[:D].norm(); y = zt[D:] / zt[D:].norm()
            Wy = pmap(torch.outer(y, y), n); Wx = pmap(torch.outer(x, x), n); lam = x @ Wy @ x
            r = ((Wy @ x - lam * x)**2).sum() + ((Wx @ y - lam * y)**2).sum(); r.backward()
            return r.item(), zt.grad.numpy()
        r = minimize(fr, rng.standard_normal(2 * D), jac=True, method='L-BFGS-B', options={'maxiter': 60000, 'gtol': 1e-16, 'ftol': 1e-24})
        if r.fun < 1e-20:
            pts.append((r.x[:D] / np.linalg.norm(r.x[:D]), r.x[D:] / np.linalg.norm(r.x[D:])))
    return pts
def thickness(v, n):
    V = v.reshape([3] * n)
    return min(np.linalg.eigvalsh(np.moveaxis(V, s, 0).reshape(3, -1) @ np.moveaxis(V, s, 0).reshape(3, -1).T)[0] for s in range(n))
def at_copy(M, s, v, n):
    V = np.moveaxis(v.reshape([3] * n), s, 0); V = np.tensordot(M, V, axes=(1, 0)); return np.moveaxis(V, 0, s).reshape(-1)
tl = []
for (i, j) in [(0, 1), (0, 2), (1, 2)]:
    E = np.zeros((3, 3)); E[i, j] = E[j, i] = 1 / np.sqrt(2); tl.append(E)
tl += [np.diag([1., -1., 0.]) / np.sqrt(2), np.diag([1., 1., -2.]) / np.sqrt(6)]
def quantities(pmap, a, n, x, y, s):
    Wy = pmap(torch.outer(torch.tensor(y), torch.tensor(y)), n).numpy(); lam = x @ Wy @ x
    Xs = np.moveaxis(x.reshape([3] * n), s, 0).reshape(3, -1); Ys = np.moveaxis(y.reshape([3] * n), s, 0).reshape(3, -1)
    tau = np.trace((Xs.T @ Xs) @ pmap(torch.tensor(Ys.T @ Ys), n - 1).numpy()); trr = np.trace((Xs @ Xs.T) @ (Ys @ Ys.T))
    return lam, tau, trr
def hess_trace(pmap, a, n, x, y, s):
    def f(u, v): return v @ pmap(torch.outer(torch.tensor(u), torch.tensor(u)), n).numpy() @ v
    def bil(p, q, r_, t_):   # <p (x) q | C | r (x) t> with C the Choi matrix of P^{(x)n}: = <q| P(p r^T sym) |t>
        M = pmap(torch.tensor(np.outer(p, r_) + np.outer(r_, p)) / 2, n).numpy(); return q @ M @ t_
    lam = f(x, y); tr = 0.0; fd = 0.0
    for E in tl:
        Hx, Hy = at_copy(E, s, x, n), at_copy(E, s, y, n)
        tr += (2 * f(Hx, y) + 2 * f(x, Hy) + 4 * bil(Hx, y, x, Hy) + 4 * bil(Hx, Hy, x, y)
               - 2 * lam * (Hx @ Hx + Hy @ Hy) - 8 * lam * (x @ Hx) * (y @ Hy))
        h = 1e-4; vals = []
        for t in (-h, 0.0, h):
            G = np.eye(3) + t * E + t * t * (E @ E) / 2 + t**3 * (E @ E @ E) / 6
            gx, gy = at_copy(G, s, x, n), at_copy(G, s, y, n)
            vals.append(f(gx, gy) / ((gx @ gx) * (gy @ gy)))
        fd += (vals[0] - 2 * vals[1] + vals[2]) / h**2
    return tr, fd
print("Part 1: trace identity")
rng = np.random.default_rng(726)
for a in (0.75, 0.7):
    pmap = make(a); n = 3
    p = np.array([1, 1j, 0]) / np.sqrt(2); g = np.kron(np.kron(p, p), p); ghz = ((g + g.conj()) / np.sqrt(2)).real
    e1 = np.zeros(3); e1[0] = 1; prod = np.kron(np.kron(e1, e1), e1)
    tests = [('aligned product', prod, prod), ('plane GHZ', ghz, ghz)] + [('thick saddle', x, y) for x, y in critical(pmap, n, rng, 3) if min(thickness(x, n), thickness(y, n)) > 1e-3][:2]
    for name, x, y in tests:
        lam, tau, trr = quantities(pmap, a, n, x, y, 0)
        pred = -4 * (lam * (4 * a - 1 + 2 * trr) - (a - .5) * (1 + 4 * a) * tau)
        tr, fd = hess_trace(pmap, a, n, x, y, 0)
        good = abs(tr - pred) < 1e-9 and abs(fd - pred) < 1e-4
        ok &= good
        print(f"  a = {a}, {name}: lam {lam:+.6f}; trace from the Hessian formula {tr:+.8f}, finite differences {fd:+.6f}, "
              f"identity {pred:+.8f}{'' if good else '  FAIL'}", flush=True)
print("Part 2: ratios at thick critical points, a = 3/4")
a = 0.75; pmap = make(a)
for n, runs in ((3, 15), (4, 8)):
    allr = []; cnt = 0
    for x, y in critical(pmap, n, rng, runs):
        if min(thickness(x, n), thickness(y, n)) < 1e-3: continue
        cnt += 1; rs = []
        for s in range(n):
            lam, tau, trr = quantities(pmap, a, n, x, y, s)
            rs.append((a - .5) * (1 + 4 * a) * tau / (lam * (4 * a - 1 + 2 * trr)))
        allr += rs
        print(f"  n = {n}: thick critical point lam = {lam:+.6f}, thickness {min(thickness(x, n), thickness(y, n)):.3f}; "
              f"ratios per copy {np.round(rs, 3)}", flush=True)
    print(f"  n = {n}: {cnt} thick critical points, largest ratio over all copies {max(allr):.4f}" if allr else f"  n = {n}: none found")
    ok &= (not allr) or max(allr) < 1
print("PASS" if ok else "FAIL")
