"""Opposite-sign route.  For orthonormal u, v in (C^3)^(x)n:
   g(u,v) = <u (x) v| (x)_i (2 F_i - I) |u (x) v> = sum_T (-1)^{n-|T|} 2^{|T|} Tr(rho^u_T rho^v_T),
   q_n(mu P_u - nu P_v) = 2^-n [mu^2 g(u,u) + nu^2 g(v,v) - 2 mu nu g(u,v)],  g(w,w) >= 1 always.
OS_n:  g(u,v) <= 1 for all orthonormal u, v   ==>  opposite-sign Hermitian endpoint at n copies.
Maximise g(u,v) by L-BFGS from random and structured starts."""
import sys, itertools, numpy as np
from scipy.optimize import minimize
n = int(sys.argv[1]); NS = int(sys.argv[2]); seed = int(sys.argv[3]); d = 3
rng = np.random.default_rng(seed); D = d ** n
subsets = [T for k in range(1, n) for T in itertools.combinations(range(n), k)]
coef = {T: (-1) ** (n - len(T)) * 2 ** len(T) for T in subsets}
c0 = (-1) ** n   # T = empty: t = 1 ; T = [n]: t = |<u|v>|^2 = 0
def red(x, T):
    X = x.reshape((d,) * n); rest = [i for i in range(n) if i not in T]
    X = np.transpose(X, list(T) + rest).reshape(d ** len(T), d ** len(rest)); return X @ X.conj().T
def gval(u, v):
    return c0 + sum(cf * np.real(np.trace(red(u, T) @ red(v, T))) for T, cf in coef.items())
def unpack(z):
    a = z[:D] + 1j * z[D:2 * D]; b = z[2 * D:3 * D] + 1j * z[3 * D:]
    u = a / np.linalg.norm(a); b = b - np.vdot(u, b) * u; v = b / np.linalg.norm(b); return u, v
def f(z):
    u, v = unpack(z); return -gval(u, v)
best = -1e9
for s in range(NS):
    z0 = rng.standard_normal(4 * D) * (1.0 if s % 2 else 0.3)
    r = minimize(f, z0, method='L-BFGS-B', options={'maxiter': 5000})
    u, v = unpack(r.x); gv = gval(u, v)
    best = max(best, gv)
    print(f"n={n} start {s}: max g(u,v) = {gv:.10f}   g(u,u) = {gval(u,u):.4f}  g(v,v) = {gval(v,v):.4f}", flush=True)
print(f"n={n}: MAX g(u,v) over {NS} starts = {best:.10f}  (OS_n claims <= 1)", flush=True)
