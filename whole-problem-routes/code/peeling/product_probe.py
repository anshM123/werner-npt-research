"""Numerical check of Proposition 7.30 (floating point).  For states rho on n qutrit copies (random mixed, random pure,
near products, balanced) and real unit product vectors v (random, and locally optimised to minimise the value), check

    <v|P_{3/4}^{(x)n}(rho)|v>  >=  (3/4)^n 3^{-sum_s <v_s|rho_s|v_s>}  >=  4^{-n} 3^{2 sum_s tau_s}  >=  prod_s (1/4 + tau_s/2).

Ends with PASS or FAIL."""
import numpy as np
from scipy.optimize import minimize

rng = np.random.default_rng(7)


def pmap(X, m):
    T = X.reshape([3] * (2 * m))
    for s in range(m):
        T = np.moveaxis(T, [s, m + s], [0, 1])
        sh = T.shape
        M = T.reshape(3, 3, -1)
        tr = np.einsum('iir->r', M)
        T = np.moveaxis((-0.25 * (M + M.transpose(1, 0, 2)) + 0.75 * np.eye(3)[:, :, None] * tr).reshape(sh), [0, 1], [s, m + s])
    return T.reshape(3 ** m, 3 ** m)


def marg(rho, m, s):
    T = np.moveaxis(rho.reshape([3] * (2 * m)), [s, m + s], [0, m]).reshape(3, 3 ** (m - 1), 3, 3 ** (m - 1))
    return np.einsum('iaja->ij', T)


def prod_vec(vs):
    out = np.array([1.0])
    for v in vs:
        out = np.kron(out, v / np.linalg.norm(v))
    return out


def states(m):
    d = 3 ** m
    out = []
    for r in (1, 2, d):
        G = rng.standard_normal((d, r))
        out.append(('random rank %d' % r, G @ G.T / np.trace(G @ G.T)))
    e = prod_vec([rng.standard_normal(3) for _ in range(m)])
    G = np.outer(e, np.ones(2)) + 0.1 * rng.standard_normal((d, 2))
    out.append(('near product', G @ G.T / np.trace(G @ G.T)))
    if m >= 2:
        om = np.eye(3).reshape(-1) / np.sqrt(3)
        y = om
        for _ in range(m - 2):
            y = np.kron(y, prod_vec([rng.standard_normal(3)]))
        out.append(('max entangled pair (x) product', np.outer(y, y)))
    return out


ok = True
worst = np.inf
for m in range(1, 6):
    for name, rho in states(m):
        W = pmap(rho, m)
        rs = [marg(rho, m, s) for s in range(m)]
        taus = [np.linalg.eigvalsh(r)[0] for r in rs]
        tc = np.prod([0.25 + 0.5 * t for t in taus])
        exp_bound = 4.0 ** -m * 3.0 ** (2 * sum(taus))

        def val(z):
            vs = [z[3 * s:3 * s + 3] for s in range(m)]
            v = prod_vec(vs)
            return v @ W @ v

        for trial in range(6):
            z0 = rng.standard_normal(3 * m)
            if trial == 0:
                z0 = np.concatenate([np.linalg.eigh(r)[1][:, -1] for r in rs])
            res = minimize(val, z0, method='BFGS')
            vs = [res.x[3 * s:3 * s + 3] / np.linalg.norm(res.x[3 * s:3 * s + 3]) for s in range(m)]
            v = prod_vec(vs)
            lhs = v @ W @ v
            jensen = 0.75 ** m * 3.0 ** (-sum(vs[s] @ rs[s] @ vs[s] for s in range(m)))
            good = lhs >= jensen - 1e-12 and jensen >= exp_bound - 1e-12 and exp_bound >= tc - 1e-12
            ok &= bool(good)
            worst = min(worst, lhs / tc)
        print(f'n = {m}, {name:32s}: min over product probes of value / TC bound = {worst:.6f}', flush=True)
        worst = np.inf
print('PASS' if ok else 'FAIL')
