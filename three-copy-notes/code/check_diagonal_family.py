"""Theorem (diagonal planes). Let V be a plane inside span{|jjj> : j = 1..d} in (C^d)^{(x)3}, and let g be the diagonal
of the rank-two projector G on C^d with P_V = sum_{jk} G_jk |jjj><kkk|. Then
    q_3(P_V) = 3/2 - (3/4) sum_j g_j^2,   lambda_3(V) = third largest g_j,
and q_3(P_V) >= (3/4) lambda_3(V), i.e. sum_j g_j^2 + g_(3) <= 2 on the Schur-Horn polytope {0 <= g <= 1, sum g = 2}.
This script checks the closed form on random planes (d = 3, 4) and the polytope inequality by sampling (d = 3..8).
Run: python check_diagonal_family.py"""
import itertools
import numpy as np

rng = np.random.default_rng(7)


def q3(X, d):
    T = X.reshape([d] * 6)
    tot = 0.0
    for J in itertools.product([0, 1], repeat=3):
        A, m = T, 3
        for i in reversed(range(3)):
            if J[i]:
                A = np.trace(A, axis1=i, axis2=i + m)
                m -= 1
        M = A.reshape(d ** m, d ** m)
        tot += (-0.5) ** sum(J) * np.vdot(M, M).real
    return tot


err, worst_margin = 0.0, np.inf
for d in (3, 4):
    for _ in range(100):
        Q, _ = np.linalg.qr(rng.standard_normal((d, 2)) + 1j * rng.standard_normal((d, 2)))
        u = np.zeros(d ** 3, complex); v = np.zeros(d ** 3, complex)
        for j in range(d):
            u[j * (d * d + d + 1)] = Q[j, 0]; v[j * (d * d + d + 1)] = Q[j, 1]
        P = np.outer(u, u.conj()) + np.outer(v, v.conj())
        g = np.real(np.diag(Q @ Q.conj().T))
        val = q3(P, d)
        err = max(err, abs(val - (1.5 - 0.75 * (g ** 2).sum())))
        rho = np.trace(P.reshape(d, d * d, d, d * d), axis1=1, axis2=3)
        lam3 = np.sort(np.linalg.eigvalsh(rho))[-3]
        worst_margin = min(worst_margin, val - 0.75 * lam3)
print(f"closed form q_3 = 3/2 - (3/4) sum g^2: max error {err:.1e};  min of q_3 - (3/4) lambda_3: {worst_margin:+.2e}")
worst = -np.inf
for d in range(3, 9):
    for _ in range(100000):
        g = 2 * rng.dirichlet(np.ones(d) * rng.choice([0.2, 0.5, 1.0, 3.0]))
        if g.max() > 1:
            continue
        s = np.sort(g)[::-1]
        worst = max(worst, (g ** 2).sum() + s[2] - 2)
print(f"max of sum g^2 + g_(3) - 2 over sampled points of the Schur-Horn polytope, d = 3..8: {worst:+.2e} (<= 0)")
