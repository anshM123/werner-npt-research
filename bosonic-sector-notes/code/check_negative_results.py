"""Section 3, "what does not work" (t = 1/2).
   (a) U(1) charge sectors (blocks E_l C E_m with l - m = delta fixed) are indefinite, n <= 4;
   (b) spaces of degree-n forms through m generic points (V_S = span{p^n : p in S}^perp) are not positive unless
       2-dimensional, n = 2, 3;
   (c) pairs gamma_k p_k + gamma_{k-1} p_{k-1} on rank-2 C: the top pair stays >= 0, lower pairs go negative, n <= 4."""
import numpy as np
from math import factorial, sqrt
from scipy.optimize import minimize
from bosonic import basis, q_and_grad, A_op, A_adj, gamma

rng = np.random.default_rng(6)

def qform(n):
    D = len(basis(n)); Q = np.zeros((D * D, D * D))
    for e in range(D * D):
        X = np.zeros(D * D); X[e] = 1
        Q[:, e] = q_and_grad(X.reshape(D, D), n)[1].ravel()
    return (Q + Q.T) / 2, D

print("(a) charge sectors: minimum eigenvalue for each delta = l - m")
for n in range(1, 5):
    Q, D = qform(n); occ = np.array([b[2] for b in basis(n)])
    out = []
    for dl in range(-n, n + 1):
        idx = [r * D + c for r in range(D) for c in range(D) if occ[r] - occ[c] == dl]
        out.append(round(float(np.linalg.eigvalsh(Q[np.ix_(idx, idx)])[0]), 5))
    print(f"  n={n}: delta = -n..n -> {out}")

print("(b) forms through m generic points: minimum eigenvalue over 30 random point sets")
for n in (2, 3):
    Q, D = qform(n); B = basis(n)
    def coherent(p):
        return np.array([sqrt(factorial(n) / (factorial(a) * factorial(b) * factorial(c))) * p[0] ** a * p[1] ** b * p[2] ** c
                         for (a, b, c) in B])
    for m in range(1, D - 1):
        worst = np.inf
        for _ in range(30):
            P = rng.standard_normal((m, 3)) + 1j * rng.standard_normal((m, 3))
            Cm = np.array([coherent(p / np.linalg.norm(p)) for p in P]).conj()
            Vb = np.linalg.svd(Cm)[2][m:].conj().T; k = Vb.shape[1]
            L = np.zeros((D * D, k * D), complex)
            for e in range(k * D):
                X = np.zeros(k * D, complex); X[e] = 1; L[:, e] = (Vb @ X.reshape(k, D)).ravel()
            M = L.conj().T @ Q @ L
            worst = min(worst, np.linalg.eigvalsh((M + M.conj().T) / 2)[0])
        print(f"  n={n} m={m} (dim {D - m}): {worst:+.4f}")

print("(c) pairs over rank-2 C (local minimisation)")
for n in (2, 3, 4):
    D = len(basis(n)); M = np.zeros((D * D, D * D))
    for e in range(D * D):
        X = np.zeros(D * D); X[e] = 1; M[:, e] = A_adj(A_op(X.reshape(D, D), n), n).ravel()
    w, V = np.linalg.eigh((M + M.T) / 2)
    P = {k: (lambda S: S @ S.T)(V[:, np.abs(w - (n - k) * (n + k + 2)) < 1e-6]) for k in range(n + 1)}
    res = []
    for k in range(n, 0, -2):
        Mop = float(gamma(n, k)) * P[k] + float(gamma(n, k - 1)) * P[k - 1]
        def f(x):
            Y = (x[:2 * D] + 1j * x[2 * D:4 * D]).reshape(D, 2); Z = (x[4 * D:6 * D] + 1j * x[6 * D:]).reshape(D, 2)
            C = Y @ Z.conj().T; nn = np.vdot(C, C).real; G = (Mop @ C.ravel()).reshape(D, D)
            val = np.vdot(C, G).real / nn; Gc = (G - val * C) / nn; GY = Gc @ Z; GZ = Gc.conj().T @ Y
            return val, 2 * np.concatenate([GY.real.ravel(), GY.imag.ravel(), GZ.real.ravel(), GZ.imag.ravel()])
        best = min(minimize(f, rng.standard_normal(8 * D), jac=True, method='L-BFGS-B',
                            options={'maxiter': 20000, 'gtol': 1e-12}).fun for _ in range(12))
        res.append((k, k - 1, round(best, 5)))
    print(f"  n={n}: (k, k-1, min) {res}")
