"""Independent numerical check of the rank-preserving realification identity (II.4):
q_{n+2r}(C_r) = 1/2[(11/12)^r q_n(C) + (-1/12)^r Re b_n(C, conj C) + (2/3)^r q_n(T)], C = U S V^dag, T = U S V^T,
C_r = sum_j s_j (R_r u_j)(R_r v_j)^T, R_r z = (z (x) F_r + conj(z) (x) G_r)/sqrt2, F_r = f^{(x)r}, G_r = conj F_r,
f = (|00> + w|11> + w^2|22>)/sqrt3; plus: C_r is real, rank 2, same norm."""
import numpy as np
rng = np.random.default_rng(7)
om = np.exp(2j * np.pi / 3)


def qform(X, Y, n):
    """b_n(X, Y) = <X, M^{(x)n}(Y)>, M(Z) = Z - Tr(Z) I/2 on each copy (X, Y operators on (C^3)^{(x)n})"""
    T = Y.reshape((3,) * (2 * n))
    for i in range(n):
        tr = np.trace(T, axis1=i, axis2=n + i)                      # partial trace over copy i
        T = T - 0.5 * np.moveaxis(np.multiply.outer(tr, np.eye(3)), [-2, -1], [i, n + i])
    return np.vdot(X.reshape(-1), T.reshape(-1))


f = np.zeros(9, complex)
for k in range(3):
    f[3 * k + k] = om ** k / np.sqrt(3)
g = f.conj()
worst = 0.0
for n, r in [(1, 1), (1, 2), (2, 1)]:
    for trial in range(3):
        D = 3 ** n
        U, _ = np.linalg.qr(rng.normal(size=(D, 2)) + 1j * rng.normal(size=(D, 2)))
        V, _ = np.linalg.qr(rng.normal(size=(D, 2)) + 1j * rng.normal(size=(D, 2)))
        s = rng.uniform(0.3, 2, 2)
        C = U @ np.diag(s) @ V.conj().T
        T = U @ np.diag(s) @ V.T
        F = f; G = g
        for _ in range(r - 1):
            F = np.kron(F, f); G = np.kron(G, g)
        R = lambda z: (np.kron(z, F) + np.kron(z.conj(), G)) / np.sqrt(2)
        Cr = sum(s[j] * np.outer(R(U[:, j]), R(V[:, j])) for j in range(2))
        # the flag copies are placed after the n original copies (copy order: n originals, then 2r flag copies)
        lhs = qform(Cr, Cr, n + 2 * r).real
        rhs = 0.5 * ((11 / 12) ** r * qform(C, C, n).real + (-1 / 12) ** r * qform(C, C.conj(), n).real
                     + (2 / 3) ** r * qform(T, T, n).real)
        real = np.abs(Cr.imag).max(); rank = np.linalg.matrix_rank(Cr, tol=1e-9)
        nrm = abs(np.linalg.norm(Cr) - np.linalg.norm(C))
        worst = max(worst, abs(lhs - rhs))
        print(f'n={n} r={r}: lhs {lhs:+.12f} rhs {rhs:+.12f} |diff| {abs(lhs-rhs):.1e}; max |Im C_r| {real:.1e}; rank {rank}; '
              f'norm diff {nrm:.1e}')
print('PASS' if worst < 1e-10 else 'FAIL')
