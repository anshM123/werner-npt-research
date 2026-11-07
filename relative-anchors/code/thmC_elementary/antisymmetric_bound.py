"""Floating-point check of the antisymmetric branch (NOTES.md, Section 3).  With free squared mass 1, restrict the form of
Theorem C to antisymmetric X (q decouples), split X into hard Z (Zw = 0) and soft w_i w_j (d_i - d_j), eliminate the
soft part, and compare lambda_min of the hard Schur operator (Frobenius metric) with h(t) - delta(t) and with 4/9,
h = 1/2 - 1/((3K+1)(6K+5)), delta = t/(8 K^2 (K+1)^3), K = t + 1/2, t = sigma.  No claim rests on it."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from freeblock import free_block
rng = np.random.default_rng(5)
worst = worst2 = 9.0
for it in range(200):
    u = int(rng.integers(3, 8))
    xU = np.exp(rng.normal(0, 1.5, u)); xU = xU / xU.sum()
    t = float(np.exp(rng.uniform(-6, 6)))
    w = np.sqrt(xU)
    C = np.array(free_block([np.sqrt(t)], list(w), float))
    Cxx = C[:-1, :-1]
    basis = []
    for i in range(u):
        for j in range(i + 1, u):
            B = np.zeros((u, u)); B[i, j] = 1; B[j, i] = -1
            basis.append(B.reshape(-1) / np.sqrt(2))
    Ba = np.array(basis).T
    soft = []
    for k in range(u):
        Mk = np.array([[w[i] * w[j] * ((i == k) - (j == k)) for j in range(u)] for i in range(u)])
        soft.append(Mk.reshape(-1))
    Qs, sv, _ = np.linalg.svd(np.array(soft).T, full_matrices=False)
    Qs = Qs[:, sv > 1e-12]
    P = Ba - Qs @ (Qs.T @ Ba)
    Uh, sv, _ = np.linalg.svd(P, full_matrices=False)
    Hb = Uh[:, sv > 1e-9]
    assert Hb.shape[1] == u * (u - 1) // 2 - (u - 1) and Qs.shape[1] == u - 1
    Fhh, Fhs, Fss = Hb.T @ Cxx @ Hb, Hb.T @ Cxx @ Qs, Qs.T @ Cxx @ Qs
    Hc = Fhh - Fhs @ np.linalg.solve(Fss, Fhs.T)
    lam = np.linalg.eigvalsh((Hc + Hc.T) / 2)[0]
    K = t + 0.5
    bound = 0.5 - 1 / ((3 * K + 1) * (6 * K + 5)) - t / (8 * K ** 2 * (K + 1) ** 3)
    worst = min(worst, lam - 4 / 9); worst2 = min(worst2, lam - bound)
print(f'200 configurations (u = 3..7, t from e^-6 to e^6): min (lambda_min(H_c,-) - 4/9) = {worst:.6f}; '
      f'min (lambda_min(H_c,-) - (h - delta)) = {worst2:.2e}')
print('PASS' if worst > 0 and worst2 > -1e-10 else 'FAIL')
