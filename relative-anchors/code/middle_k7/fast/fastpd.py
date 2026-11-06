"""Exact positive-definiteness test by recursive Schur complements in python-flint (fmpq_mat), with the same result as
the sequential symmetric elimination: (True, None) if every pivot is > 0, else (False, index of the first pivot <= 0).
The pivots of M in natural order are those of its leading block A followed by those of the Schur complement
C - B^T A^{-1} B, so the recursion visits the same pivots."""
from fractions import Fraction as Fr
import flint


def _sub(ent, n, r0, r1, c0, c1):
    return flint.fmpq_mat(r1 - r0, c1 - c0, [ent[i * n + j] for i in range(r0, r1) for j in range(c0, c1)])


def _rec(A, off):
    n = A.nrows()
    ent = A.entries()
    if n <= 24:
        M = [[ent[i * n + j] for j in range(n)] for i in range(n)]
        for c in range(n):
            p = M[c][c]
            if p <= 0:
                return False, off + c
            for i in range(c + 1, n):
                if M[i][c] != 0:
                    f = M[i][c] / p
                    Mi, Mc = M[i], M[c]
                    for l in range(c + 1, n):
                        if Mc[l] != 0:
                            Mi[l] -= f * Mc[l]
        return True, None
    h = n // 2
    A11 = _sub(ent, n, 0, h, 0, h)
    ok, c = _rec(A11, off)
    if not ok:
        return ok, c
    A12 = _sub(ent, n, 0, h, h, n)
    A22 = _sub(ent, n, h, n, h, n)
    S = A22 - A12.transpose() * A11.solve(A12)
    return _rec(S, off + h)


def ldl_pd_fast(M):
    n = len(M)
    if n == 0:
        return True, None
    A = flint.fmpq_mat(n, n, [flint.fmpq(Fr(x).numerator, Fr(x).denominator) for row in M for x in row])
    return _rec(A, 0)
