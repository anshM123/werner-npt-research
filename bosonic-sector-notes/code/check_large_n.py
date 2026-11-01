"""Bosonic sector at large n: minimise q_n(C)/||C||^2 over rank-2 C = Y Z^* in End(Sym^n C^3).
q_n(C) = sum_j c_j ||A^j C||^2, c_j = (-1/2)^j (n-j)!/(j! n!), A(C) = sum_i a_i C a_i^dag (bosonic annihilators),
the formula of Section 1 (checked against direct partial traces by check_large_n_formula.py).  Sparse annihilators.
Starts: random; near qubit zeros (C supported on Sym^n(span{e0,e1}) with SU(2) spins l != n mod 2, plus leakage into
level 2); Dicke / GHZ-type basis pairs.  usage: python check_large_n.py n starts seed   (the runs in logs/check_large_n.log used seed = n, 18 starts)"""
import sys, numpy as np, scipy.sparse as sp
from math import lgamma, exp
from scipy.optimize import minimize
n, NS, SEED = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
rng = np.random.default_rng(SEED)
def basis(m): return [(a, b, m - a - b) for a in range(m + 1) for b in range(m + 1 - a)]
ANN = {}
def ann(m):
    if m in ANN: return ANN[m]
    B, Bm = basis(m), basis(m - 1); idx = {t: k for k, t in enumerate(Bm)}; out = []
    for i in range(3):
        r, c, v = [], [], []
        for k, t in enumerate(B):
            if t[i] > 0:
                tt = list(t); tt[i] -= 1; r.append(idx[tuple(tt)]); c.append(k); v.append(np.sqrt(t[i]))
        A = sp.csr_matrix((v, (r, c)), shape=(len(Bm), len(B)))
        out.append((A, A.T.tocsr()))
    ANN[m] = out; return out
def Aop(X, m): return sum(a @ (a @ X.T).T for a, aT in ann(m))      # a X a^T  (a real)
def Adag(X, m): return sum(aT @ (aT @ X.T).T for a, aT in ann(m))   # a^T X a
def cj(j): return (-0.5) ** j * exp(lgamma(n - j + 1) - lgamma(j + 1) - lgamma(n + 1))
CJ = [cj(j) for j in range(n + 1)]
def q_and_grad(C):
    Xs = [C]
    for j in range(1, n + 1): Xs.append(Aop(Xs[-1], n - j + 1))
    q = sum(CJ[j] * np.vdot(X, X).real for j, X in enumerate(Xs))
    G = CJ[n] * Xs[n]
    for j in range(n, 0, -1): G = CJ[j - 1] * Xs[j - 1] + Adag(G, n - j + 1)
    return q, G
B = basis(n); D = len(B); pos = {t: k for k, t in enumerate(B)}
def f(x):
    Y = (x[:2 * D] + 1j * x[2 * D:4 * D]).reshape(D, 2); Z = (x[4 * D:6 * D] + 1j * x[6 * D:]).reshape(D, 2)
    C = Y @ Z.conj().T; nn = np.vdot(C, C).real; q, G = q_and_grad(C)
    val = q / nn; Gc = (G - val * C) / nn
    GY = Gc @ Z; GZ = Gc.conj().T @ Y
    return val, 2 * np.concatenate([GY.real.ravel(), GY.imag.ravel(), GZ.real.ravel(), GZ.imag.ravel()])
def pack(Y, Z): return np.concatenate([Y.real.ravel(), Y.imag.ravel(), Z.real.ravel(), Z.imag.ravel()])
qubit = [k for k, t in enumerate(B) if t[2] == 0]           # Sym^n(span{e0, e1})
def start(kind):
    if kind == 'rand':
        return rng.standard_normal(8 * D)
    Y = np.zeros((D, 2), complex); Z = np.zeros((D, 2), complex)
    if kind == 'qubit':
        Y[qubit] = rng.standard_normal((len(qubit), 2)) + 1j * rng.standard_normal((len(qubit), 2))
        Z[qubit] = rng.standard_normal((len(qubit), 2)) + 1j * rng.standard_normal((len(qubit), 2))
        eps = 10 ** rng.uniform(-3, -0.5)
        Y += eps * (rng.standard_normal((D, 2)) + 1j * rng.standard_normal((D, 2)))
        Z += eps * (rng.standard_normal((D, 2)) + 1j * rng.standard_normal((D, 2)))
    else:  # 'ghz': three-level GHZ / Dicke pairs plus noise
        ks = rng.choice(D, size=4, replace=False)
        for r in range(2):
            Y[ks[r], r] = 1; Z[ks[2 + r] if rng.random() < 0.5 else ks[r], r] = 1
        Y[pos[(n, 0, 0)], 0] += 1; Y[pos[(0, n, 0)], 1] += 1; Z[pos[(n, 0, 0)], 0] += 1; Z[pos[(0, 0, n)], 1] += 1
        Y += 0.05 * (rng.standard_normal((D, 2)) + 1j * rng.standard_normal((D, 2)))
        Z += 0.05 * (rng.standard_normal((D, 2)) + 1j * rng.standard_normal((D, 2)))
    return pack(Y, Z)
kinds = (['rand', 'qubit', 'ghz'] * NS)[:NS]
best = (np.inf, None)
for s, kind in enumerate(kinds):
    r = minimize(f, start(kind), jac=True, method='L-BFGS-B', options={'maxiter': 6000, 'gtol': 1e-11, 'ftol': 1e-15})
    if r.fun < best[0]: best = (r.fun, kind)
    print(f"n={n} D={D} start {s:2d} ({kind:5s}): q/||C||^2 = {r.fun:+.3e}  its {r.nit}", flush=True)
print(f"n={n}: MIN over {NS} starts = {best[0]:+.3e} ({best[1]})", flush=True)
