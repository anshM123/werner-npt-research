"""Independent check of the cyclic qutrit code result (II.12) from the definition q_n(C) = <vec C, M^{(x)n} vec C>,
M = I - |Omega><Omega|/2 on C^3_row (x) C^3_col: least eigenvalue of the compressed form for C with range in Q_n
(one-sided, arbitrary co-range) and for C = P C P (two-sided), n = 4 (dense) and n = 5 (Lanczos).  Floating point.
Expected (exact, from the transfer formula): n = 4: 144/432; n = 5: 639/2592 (one-sided); two-sided n = 4: 144/432."""
import itertools, sys
import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh
om = np.exp(2j * np.pi / 3)
Mloc = np.eye(9).reshape(3, 3, 3, 3) - 0.5 * np.einsum('ab,cd->abcd', np.eye(3), np.eye(3))   # M[(a,b),(a',b')]


def code_vectors(n):
    xs = np.array(list(itertools.product(range(3), repeat=n)))
    ph = sum(xs[:, i] * xs[:, (i + 1) % n] for i in range(n))
    tot = xs.sum(axis=1)
    return [om ** ((ph + t * tot) % 3) / 3 ** (n / 2) for t in range(3)]     # each of length 3^n


def apply_M(vec, n):
    """vec: row (x) col vector of length 3^(2n), indices (x_1..x_n, y_1..y_n); apply M on each pair (x_i, y_i)"""
    T = vec.reshape((3,) * (2 * n))
    for i in range(n):
        T = np.tensordot(Mloc, T, axes=([2, 3], [i, n + i]))      # new axes (a, b) in front
        T = np.moveaxis(T, [0, 1], [i, n + i])
    return T.reshape(-1)


def one_sided(n):
    g = code_vectors(n)
    dimH = 3 ** n
    def W(c):          # c: 3 x dimH coefficients -> vec C = sum_t g_t (x) c_t
        c = c.reshape(3, dimH)
        return sum(np.kron(g[t], c[t]) for t in range(3))
    def Wh(v):
        V = v.reshape(dimH, dimH)
        return np.concatenate([g[t].conj() @ V for t in range(3)])
    op = LinearOperator((3 * dimH, 3 * dimH), matvec=lambda c: Wh(apply_M(W(c), n)), dtype=complex)
    if n <= 4:
        Mat = np.column_stack([op.matvec(e) for e in np.eye(3 * dimH)])
        return np.linalg.eigvalsh((Mat + Mat.conj().T) / 2)[0]
    return eigsh(op, k=1, which='SA', tol=1e-12)[0][0]


def two_sided(n):
    g = code_vectors(n)
    vecs = [np.kron(g[s], g[t].conj()) for s in range(3) for t in range(3)]
    Mat = np.array([[vecs[i].conj() @ apply_M(vecs[j], n) for j in range(9)] for i in range(9)])
    return np.linalg.eigvalsh((Mat + Mat.conj().T) / 2)[0]


print(f'n = 4: one-sided least eigenvalue {one_sided(4):.12f} (144/432 = {144/432:.12f}); two-sided {two_sided(4):.12f}')
print(f'n = 5: one-sided least eigenvalue {one_sided(5):.12f} (639/2592 = {639/2592:.12f}); two-sided {two_sided(5):.12f} '
      f'(Theorem A, n = 5: {639/2592:.12f})')
