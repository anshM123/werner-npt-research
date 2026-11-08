"""Weighted Plucker bound (audit 42): for orthonormal u, v on n qutrits,
   alpha = <uu|D|uu>, beta = <vv|D|vv>, gamma = <uv|D|uv>,  D = prod_s(2F_s - I) - 2F_all + I (two replicas),
   sharp endpoint <=> gamma >= -sqrt(alpha beta).
Cheapest uniform certificate: G = D (PSD on Sym^2): suffices that gamma + |<uu|D|vv>| >= 0.
Quantities via marginals: <uv|F_S|uv> = Tr(rho^u_S rho^v_S), <uu|F_S|uu> = Tr((rho^u_S)^2),
<uu|F_S|vv> = Tr[(Tr_{S^c}|v><u|)^2]."""
import itertools, sys
import numpy as np
from scipy.optimize import minimize
n = int(sys.argv[1]); d = 3; N = d ** n
subsets = [S for k in range(n + 1) for S in itertools.combinations(range(n), k)]
coef = {}
for S in subsets:
    c = (-1) ** (n - len(S)) * 2 ** len(S)
    if len(S) == n:
        c -= 2
    if len(S) == 0:
        c += 1
    coef[S] = c
def ptrace_keep(M, S):
    """partial trace of an N x N operator over the complement of S (keep S)"""
    T = M.reshape([d] * (2 * n))
    rest = [t for t in range(n) if t not in S]
    for k, t in enumerate(sorted(rest, reverse=True)):
        T = np.trace(T, axis1=t, axis2=t + T.ndim // 2)
    k = len(S)
    return T.reshape(d ** k, d ** k)
def quantities(u, v):
    Pu, Pv, Mvu = np.outer(u, u.conj()), np.outer(v, v.conj()), np.outer(v, u.conj())
    al = be = ga = cr = 0.0
    for S in subsets:
        ru, rv, m = ptrace_keep(Pu, S), ptrace_keep(Pv, S), ptrace_keep(Mvu, S)
        c = coef[S]
        al += c * np.real(np.trace(ru @ ru)); be += c * np.real(np.trace(rv @ rv))
        ga += c * np.real(np.trace(ru @ rv)); cr += c * np.trace(m @ m)
    return al, be, ga, cr
def orth(x):
    u = x[:N] + 1j * x[N:2 * N]; v = x[2 * N:3 * N] + 1j * x[3 * N:]
    u = u / np.linalg.norm(u); v = v - (u.conj() @ v) * u; v = v / np.linalg.norm(v)
    return u, v
def obj_sharp(x):
    al, be, ga, cr = quantities(*orth(x))
    return ga + np.sqrt(max(al, 0) * max(be, 0))
def obj_G(x):
    al, be, ga, cr = quantities(*orth(x))
    return ga + abs(cr)
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
worst_s, worst_g = np.inf, np.inf
for st in range(int(sys.argv[3]) if len(sys.argv) > 3 else 6):
    x0 = rng.standard_normal(4 * N)
    rs = minimize(obj_sharp, x0, method='Powell', options={'maxiter': 3, 'xtol': 1e-6, 'ftol': 1e-10})
    rg = minimize(obj_G, x0, method='Powell', options={'maxiter': 3, 'xtol': 1e-6, 'ftol': 1e-10})
    worst_s, worst_g = min(worst_s, rs.fun), min(worst_g, rg.fun)
    print(f'start {st}: min of gamma + sqrt(alpha beta) = {rs.fun:+.5f}; min of gamma + |<uu|D|vv>| = {rg.fun:+.5f}', flush=True)
print(f'n = {n}: sharp bound worst {worst_s:+.5f} (should be >= 0); G = D certificate worst {worst_g:+.5f}')
