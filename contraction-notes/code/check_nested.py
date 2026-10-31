"""Corollary 3: random members of the nested class K_n (aligned at one copy, tail operator in K_{n-1}, K_1 = rank two
   on one copy).  Prints the smallest q_n/||C||^2 at t = 1/2 and t = 0.4 (nonnegative by Corollary 3), and for the last
   member the ranks of the local reductions of the range and co-range projectors at every copy together with
   tau_s = Tr(P_s Q_s), P_s and Q_s the supports of those reductions (tau_s > 2 means the copy is thick)."""
import numpy as np
from werner import d, q

rng = np.random.default_rng(7)


def cv(k):
    return rng.standard_normal(d ** k) + 1j * rng.standard_normal(d ** k)


def nested(n):
    """(Y, Z) with C = Y Z^dagger in K_n; the columns are the range and co-range vectors y_k, z_k"""
    if n == 1:
        return (rng.standard_normal((d, 2)) + 1j * rng.standard_normal((d, 2)),
                rng.standard_normal((d, 2)) + 1j * rng.standard_normal((d, 2)))
    Yt, Zt = nested(n - 1)
    a, b = cv(n - 1), cv(n - 1)
    Y = np.zeros((d ** n, 2), complex); Z = np.zeros((d ** n, 2), complex)
    for k in range(2):
        Y[:, k] = np.kron(np.eye(d)[k], a) + np.kron(np.eye(d)[2], Yt[:, k])
        Z[:, k] = np.kron(np.eye(d)[k], b) + np.kron(np.eye(d)[2], Zt[:, k])
    return Y, Z


def local_data(Y, Z, n):
    out = []
    for s in range(n):
        def marg(V):
            P = V @ np.linalg.pinv(V.conj().T @ V) @ V.conj().T
            T = P.reshape([d] * (2 * n))
            for k in sorted([k for k in range(n) if k != s], reverse=True):
                T = np.trace(T, axis1=k, axis2=k + T.ndim // 2)
            return T
        Pr, Pc = marg(Y), marg(Z)
        er, Ur = np.linalg.eigh(Pr); ec, Uc = np.linalg.eigh(Pc)
        Ps, Qs = Ur[:, er > 1e-8], Uc[:, ec > 1e-8]
        out.append((Ps.shape[1], Qs.shape[1], round(float(np.linalg.norm(Ps.conj().T @ Qs) ** 2), 4)))
    return out


for n in (2, 3, 4):
    low = {0.5: np.inf, 0.4: np.inf}
    for trial in range(10):
        Y, Z = nested(n); C = Y @ Z.conj().T; nC = np.vdot(C, C).real
        for t in low:
            low[t] = min(low[t], q(C, n, t) / nC)
    print(f"n={n}: min q_n/||C||^2 over 10 members: t=1/2 {low[0.5]:+.4f}, t=0.4 {low[0.4]:+.4f};  "
          f"(rank of range reduction, rank of co-range reduction, tau) per copy: {local_data(Y, Z, n)}")
