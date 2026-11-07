"""Checks for Observations 7.26 and 7.27 (thick inputs) for the map P(X) = (3/4) Tr(X) I - (X + X^T)/4.

Part 1 (exact).  At n = 4 let p1 = Omega_01 Omega_23, p2 = Omega_02 Omega_13, p3 = Omega_03 Omega_12 (products of
unnormalised maximally entangled vectors), x = p3 - p1, y = p2 - c (p1 + p3) with 36 c^2 + 142 c - 49 = 0, c > 0.
Checked in exact arithmetic (polynomials in c reduced modulo 36 c^2 + 142 c - 49): P^{(x)4}(yy^T) x = lam |y|^2 x and
P^{(x)4}(xx^T) y = lam |x|^2 y with lam = (151 - sqrt 6805)/1536, and every one-copy marginal of x and of y is I/3.
Part 2 (floating point).  At that pair: lam is the least eigenvalue of P^{(x)4}(yy^T); the best one-copy squeeze margin
max_w [lam - (1/4) f_3(<w|x, <w|y)/(|<w|x|^2 |<w|y|^2)] is negative at every copy; the squeeze Hessian over all
traceless directions at all four copies is positive definite; the full Hessian on the product of spheres has negative
eigenvalues.
Part 3 (floating point).  Least value over pairs of SO(3)-invariant vectors (spin singlets) for n = 2, ..., 7.
Part 4 (floating point).  The bound f_n(x, y) >= prod_s (1/4 + kappa tau_s) (tau_s = the smaller least eigenvalue of
the one-copy marginals of x and y at copy s): adversarial searches with kappa = 1/2 at n = 2, 3 (no violation), the
singlet pair against kappa = 0.7 (violated), and the least balanced pair at n = 2 (value 11/36) against kappa = 1
(violated).
Ends with PASS or FAIL."""
import itertools
import numpy as np
import sympy as sp
import scipy.sparse as ssp
from scipy.optimize import minimize
import torch

ok = True
a = 0.75

# ------------------------------------------------------------------ part 1: exact
c = sp.symbols('c')
minpoly = sp.Poly(36 * c**2 + 142 * c - 49, c)


def vec(pairs):
    v = {}
    for idx in itertools.product(range(3), repeat=4):
        if all(idx[i] == idx[j] for i, j in pairs):
            v[idx] = sp.Integer(1)
    return v


P1, P2, P3 = vec([(0, 1), (2, 3)]), vec([(0, 2), (1, 3)]), vec([(0, 3), (1, 2)])


def lin(*terms):
    out = {}
    for coef, v in terms:
        for k, val in v.items():
            out[k] = sp.expand(out.get(k, 0) + coef * val)
    return {k: val for k, val in out.items() if val != 0}


X = lin((1, P3), (-1, P1))
Y = lin((1, P2), (-c, P1), (-c, P3))
A = sp.Rational(3, 4)


def kern(i, j, k, l):       # <k| P(|i><j|) |l>
    return A * int(i == j) * int(k == l) - sp.Rational(1, 4) * (int(i == k) * int(j == l) + int(i == l) * int(j == k))


def apply_W(v, u):          # (P^{(x)4}(v v^T) u)_m for every m, as a dict
    out = {}
    for m in itertools.product(range(3), repeat=4):
        tot = 0
        for (i, vi), (j, vj) in itertools.product(v.items(), v.items()):
            for l, ul in u.items():
                pr = 1
                for s in range(4):
                    pr *= kern(i[s], j[s], m[s], l[s])
                    if pr == 0:
                        break
                if pr != 0:
                    tot += vi * vj * ul * pr
        if tot != 0:
            out[m] = sp.expand(tot)
    return out


def reduce(e):
    return sp.Poly(sp.expand(e), c).rem(minpoly).as_expr()


def dot(u, v):
    return sp.expand(sum(u[k] * v.get(k, 0) for k in u))


nx, ny = dot(X, X), dot(Y, Y)
WyX = apply_W(Y, X)
WxY = apply_W(X, Y)
lam_num = dot(X, WyX)                     # <x| W_y |x> = lam |x|^2 |y|^2
cval = (sp.sqrt(6805) - 71) / 36
target = (151 - sp.sqrt(6805)) / 1536
lam_ok = sp.simplify(lam_num.subs(c, cval) / (nx * ny.subs(c, cval)) - target) == 0
grid4 = list(itertools.product(range(3), repeat=4))
res1 = all(reduce(WyX.get(m, 0) * nx - lam_num * X.get(m, 0)) == 0 for m in grid4)
res2 = all(reduce(WxY.get(m, 0) * ny - lam_num * Y.get(m, 0)) == 0 for m in grid4)


def marg(v, s):
    M = sp.zeros(3, 3)
    for k1, v1 in v.items():
        for k2, v2 in v.items():
            if all(k1[t] == k2[t] for t in range(4) if t != s):
                M[k1[s], k2[s]] += v1 * v2
    return M


bal = all((marg(X, s) - nx / 3 * sp.eye(3)).applyfunc(reduce) == sp.zeros(3, 3) for s in range(4)) and \
    all((marg(Y, s) - ny / 3 * sp.eye(3)).applyfunc(reduce) == sp.zeros(3, 3) for s in range(4))
print(f"Part 1 (exact): value (151 - sqrt 6805)/1536: {lam_ok}; P(yy^T)x = lam x: {res1}; P(xx^T)y = lam y: {res2}; "
      f"all one-copy marginals I/3: {bal}", flush=True)
ok &= bool(lam_ok and res1 and res2 and bal)


# ------------------------------------------------------------------ part 2: floating point at the pair
def pmap(Xm, m):
    T = Xm.reshape([3] * (2 * m))
    for s in range(m):
        T = np.moveaxis(T, [s, m + s], [0, 1])
        sh = T.shape
        Mi = T.reshape(3, 3, -1)
        tr = np.einsum('iir->r', Mi)
        out = -0.25 * (Mi + Mi.transpose(1, 0, 2))
        out += a * np.eye(3)[:, :, None] * tr
        T = np.moveaxis(out.reshape(sh), [0, 1], [s, m + s])
    return T.reshape(3**m, 3**m)


cf = float(cval)
lamf = float(target)


def dense(v):
    out = np.zeros(81)
    for k, val in v.items():
        out[np.ravel_multi_index(k, (3, 3, 3, 3))] = float(sp.sympify(val).subs(c, cf))
    return out / np.linalg.norm(out)


x, y = dense(X), dense(Y)
Wy = pmap(np.outer(y, y), 4)
lmin = np.linalg.eigvalsh((Wy + Wy.T) / 2)[0]
margins = []
for s in range(4):
    Xs = np.moveaxis(x.reshape(3, 3, 3, 3), s, 0).reshape(3, -1)
    Ys = np.moveaxis(y.reshape(3, 3, 3, 3), s, 0).reshape(3, -1)

    def neg(t):
        w = np.array([np.sin(t[0]) * np.cos(t[1]), np.sin(t[0]) * np.sin(t[1]), np.cos(t[0])])
        xi, eta = w @ Xs, w @ Ys
        return -(lamf - 0.25 * (eta @ pmap(np.outer(xi, xi), 3) @ eta) / ((xi @ xi) * (eta @ eta)))
    best = -np.inf
    for t0 in [(th, ph) for th in np.linspace(0.15, 1.45, 5) for ph in np.linspace(0, 6.0, 8)]:
        best = max(best, -minimize(neg, t0, method='Nelder-Mead', options={'xatol': 1e-10, 'fatol': 1e-14, 'maxiter': 300}).fun)
    margins.append(best)
D = 81
Wx = pmap(np.outer(x, x), 4)
Hxx = Wy + Wy.T
Hyy = Wx + Wx.T
Hxy = np.zeros((D, D))
for b in range(D):
    e = np.zeros(D)
    e[b] = 1
    M = pmap(np.outer(x, e) + np.outer(e, x), 4)
    Hxy[:, b] = (M + M.T) @ y
H = np.block([[Hxx, Hxy.T], [Hxy, Hyy]])       # Hxy[:, b] is d/dx_b of the y-gradient: the (y, x) block
Px = np.eye(D) - np.outer(x, x)
Py = np.eye(D) - np.outer(y, y)
P = np.block([[Px, np.zeros((D, D))], [np.zeros((D, D)), Py]])
Hr = P @ H @ P - 2 * lamf * P
wfull = np.linalg.eigvalsh((Hr + Hr.T) / 2)
tl = []
for (i, j) in [(0, 1), (0, 2), (1, 2)]:
    E = np.zeros((3, 3))
    E[i, j] = E[j, i] = 1 / np.sqrt(2)
    tl.append(E)
tl += [np.diag([1., -1., 0.]) / np.sqrt(2), np.diag([1., 1., -2.]) / np.sqrt(6)]


def at_copy(M, s, v):
    V = np.moveaxis(v.reshape(3, 3, 3, 3), s, 0)
    V = np.tensordot(M, V, axes=(1, 0))
    return np.moveaxis(V, 0, s).reshape(-1)


Dm = np.array([P @ np.concatenate([at_copy(E, s, x), at_copy(E, s, y)]) for s in range(4) for E in tl])
Qsq = Dm @ Hr @ Dm.T
wsq = np.linalg.eigvalsh((Qsq + Qsq.T) / 2)
p2 = abs(lmin - lamf) < 1e-12 and max(margins) < 0 and wsq[0] > 0 and wfull[0] < -1e-6
print(f"Part 2: least eigenvalue of P(yy^T) {lmin:.12f} (lam {lamf:.12f}); best one-copy squeeze margins "
      f"{np.round(margins, 6)}; joint squeeze Hessian least eigenvalue {wsq[0]:+.6f}; full Hessian: "
      f"{np.sum(wfull < -1e-9)} negative eigenvalues, least {wfull[0]:+.6f}{'' if p2 else '  FAIL'}", flush=True)
ok &= bool(p2)


# ------------------------------------------------------------------ part 3: singlet sector
def singlet_basis(n):
    L = []
    for k in range(3):
        M = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                if len({k, i, j}) == 3:
                    M[i, j] = np.linalg.det(np.eye(3)[[k, i, j]])
        L.append(M)
    Dn = 3**n
    Cas = ssp.csr_matrix((Dn, Dn))
    for k in range(3):
        tot = ssp.csr_matrix((Dn, Dn))
        for s in range(n):
            op = ssp.identity(1, format='csr')
            for t in range(n):
                op = ssp.kron(op, ssp.csr_matrix(L[k]) if t == s else ssp.identity(3, format='csr'), format='csr')
            tot = tot + op
        Cas = Cas + tot.T @ tot
    w, V = np.linalg.eigh(Cas.toarray())
    return V[:, w < 1e-8]


print("Part 3: least value over pairs of spin singlets", flush=True)
rng = np.random.default_rng(727)
for n in range(2, 8):
    B = singlet_basis(n)
    best = np.inf
    for st in range(6):
        yv = B @ rng.standard_normal(B.shape[1])
        yv /= np.linalg.norm(yv)
        prev = None
        for it in range(300):
            M = B.T @ pmap(np.outer(yv, yv), n) @ B
            wv, V = np.linalg.eigh((M + M.T) / 2)
            xv = B @ V[:, 0]
            M2 = B.T @ pmap(np.outer(xv, xv), n) @ B
            wv2, V2 = np.linalg.eigh((M2 + M2.T) / 2)
            yv = B @ V2[:, 0]
            if prev is not None and abs(prev - wv2[0]) < 1e-15:
                break
            prev = wv2[0]
        best = min(best, wv2[0])
    print(f"  n = {n}: singlet dimension {B.shape[1]}, least value {best:.10f}, ratio to 4^-n {best * 4**n:.3f}, "
          f"per-copy factor {best**(1 / n):.4f}", flush=True)
    ok &= bool(best > 4.0**-n)

# ------------------------------------------------------------------ part 4: the thickness bonus
torch.set_default_dtype(torch.float64)


def tpmap(Xm, m):
    T = Xm.reshape([3] * (2 * m))
    for s in range(m):
        T = torch.movedim(T, (s, m + s), (0, 1))
        sh = T.shape
        Mi = T.reshape(3, 3, -1)
        tr = torch.einsum('iir->r', Mi)
        out = -0.25 * (Mi + Mi.transpose(0, 1)) + a * torch.eye(3).unsqueeze(-1) * tr
        T = torch.movedim(out.reshape(sh), (0, 1), (s, m + s))
    return T.reshape(3**m, 3**m)


def tmarg(v, n):
    rows = []
    for s in range(n):
        V = torch.movedim(v.reshape([3] * n), s, 0).reshape(3, -1)
        rows.append(torch.linalg.eigvalsh(V @ V.T)[0])
    return torch.stack(rows)


def slack_search(n, kappa, starts):
    Dn = 3**n
    worst = np.inf
    e1 = np.zeros(3)
    e1[0] = 1
    prod = e1
    for _ in range(n - 1):
        prod = np.kron(prod, e1)
    for st in range(starts):
        if st % 2:
            z0 = rng.standard_normal(2 * Dn)
        else:
            z0 = np.concatenate([prod, prod]) + 0.1 * rng.standard_normal(2 * Dn)

        def obj(z):
            zt = torch.tensor(z, requires_grad=True)
            xx = zt[:Dn] / zt[:Dn].norm()
            yy = zt[Dn:] / zt[Dn:].norm()
            fv = yy @ tpmap(torch.outer(xx, xx), n) @ yy
            tau = torch.minimum(tmarg(xx, n), tmarg(yy, n)).clamp(min=0)
            val = fv - torch.prod(0.25 + kappa * tau)
            val.backward()
            return val.item(), zt.grad.numpy()
        worst = min(worst, minimize(obj, z0, jac=True, method='L-BFGS-B', options={'maxiter': 3000}).fun)
    return worst


s2, s3 = slack_search(2, 0.5, 10), slack_search(3, 0.5, 8)
viol07 = lamf - (0.25 + 0.7 / 3)**4
def balanced_min2(starts):  # least value over pairs at n = 2 with both one-copy marginals I/3 (penalty method)
    bestv = np.inf
    for st in range(starts):
        z = rng.standard_normal(18)
        for mu in (10.0, 100.0, 1000.0, 1e4):
            def objb(zz):
                zt = torch.tensor(zz, requires_grad=True); xx = zt[:9] / zt[:9].norm(); yy = zt[9:] / zt[9:].norm()
                fv = yy @ tpmap(torch.outer(xx, xx), 2) @ yy
                d = 0.0
                for v in (xx, yy):
                    for s_ in range(2):
                        V = torch.movedim(v.reshape(3, 3), s_, 0)
                        d = d + ((V @ V.T - torch.eye(3) / 3)**2).sum()
                val = fv + mu * d; val.backward()
                return val.item(), zt.grad.numpy()
            z = minimize(objb, z, jac=True, method='L-BFGS-B', options={'maxiter': 3000}).x
        xx = z[:9] / np.linalg.norm(z[:9]); yy = z[9:] / np.linalg.norm(z[9:])
        bestv = min(bestv, yy @ pmap(np.outer(xx, xx), 2) @ yy)
    return bestv
fb = balanced_min2(6)                    # balanced pairs at n = 2 (the least value is 11/36)
viol1 = fb - (0.25 + 1 / 3)**2
p4 = s2 > -1e-9 and s3 > -1e-9 and viol07 < 0 and viol1 < 0
print(f"Part 4: kappa = 1/2: least slack found {s2:+.2e} (n = 2), {s3:+.2e} (n = 3); singlet pair against kappa = 0.7: "
      f"{viol07:+.5f}; least balanced pair at n = 2 (value {fb:.6f}) against kappa = 1: {viol1:+.5f}{'' if p4 else '  FAIL'}")
ok &= bool(p4)
print("PASS" if ok else "FAIL")
