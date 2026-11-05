"""Exact rational LP for the uniform-code enumerator relaxation (python-flint fmpq, dense two-phase simplex,
Bland's rule). Variables: A_j >= 0 (j = 0..n). Equalities: A_0 = 1, sum A = 3^n/2, sum j A_j = 4n 3^{n-2},
A_j = 0 and B_j = 0 for 1 <= j <= t. Inequalities: B_m - A_m >= 0, shadow S_m >= 0.
Objective: minimise F = sum_j (-1/2)^{n-j} A_j ; q_n = (4/3^n) F."""
import sys
from math import comb
from flint import fmpq


def kraw(n, p, q, r, s):
    M = [[0] * (n + 1) for _ in range(n + 1)]
    for j in range(n + 1):
        for h in range(n - j + 1):
            c1 = comb(n - j, h) * p ** (n - j - h) * q ** h
            for g in range(j + 1):
                M[h + g][j] += c1 * comb(j, g) * r ** (j - g) * s ** g
    return M


def simplex(A, b, c):
    """min c x s.t. A x = b, x >= 0 (b >= 0 after sign fix). Two-phase with artificials. Returns (status, value, x)."""
    m, nv = len(A), len(c)
    A = [row[:] for row in A]; b = b[:]
    for i in range(m):
        if b[i] < 0:
            A[i] = [-v for v in A[i]]; b[i] = -b[i]
    # tableau with artificials
    N = nv + m
    T = [A[i] + [fmpq(1) if k == i else fmpq(0) for k in range(m)] + [b[i]] for i in range(m)]
    basis = [nv + i for i in range(m)]

    def run(cost):
        # reduced costs row
        while True:
            z = [fmpq(0)] * (N + 1)
            for i in range(m):
                cb = cost[basis[i]]
                if cb != 0:
                    row = T[i]
                    for k in range(N + 1):
                        if row[k] != 0:
                            z[k] += cb * row[k]
            enter = -1
            for k in range(N):
                if cost[k] - z[k] < 0 and k not in basis:
                    enter = k; break                      # Bland
            if enter < 0:
                return 'optimal', z[N]
            leave, best = -1, None
            for i in range(m):
                if T[i][enter] > 0:
                    ratio = T[i][N] / T[i][enter]
                    if best is None or ratio < best or (ratio == best and basis[i] < basis[leave]):
                        best, leave = ratio, i
            if leave < 0:
                return 'unbounded', None
            piv = T[leave][enter]
            T[leave] = [v / piv for v in T[leave]]
            for i in range(m):
                if i != leave and T[i][enter] != 0:
                    f = T[i][enter]
                    Ti, Tl = T[i], T[leave]
                    T[i] = [Ti[k] - f * Tl[k] for k in range(N + 1)]
            basis[leave] = enter

    cost1 = [fmpq(0)] * nv + [fmpq(1)] * m
    st, v = run(cost1)
    if v != 0:
        return 'infeasible', None, None
    # drive artificials out of the basis where possible
    for i in range(m):
        if basis[i] >= nv:
            for k in range(nv):
                if T[i][k] != 0:
                    piv = T[i][k]
                    T[i] = [vv / piv for vv in T[i]]
                    for i2 in range(m):
                        if i2 != i and T[i2][k] != 0:
                            f = T[i2][k]
                            T[i2] = [T[i2][kk] - f * T[i][kk] for kk in range(N + 1)]
                    basis[i] = k
                    break
    cost2 = list(c) + [fmpq(10) ** 30] * m          # artificials effectively forbidden
    for i in range(m):
        pass
    # remove artificial columns from consideration by large cost
    st, v = run(cost2)
    x = [fmpq(0)] * nv
    for i in range(m):
        if basis[i] < nv:
            x[basis[i]] = T[i][N]
    return st, sum(c[k] * x[k] for k in range(nv)), x


def lp(n, t):
    K = kraw(n, 1, 8, 1, -1)           # 3^n B_m / 2 = sum_j K[m][j] A_j
    Ks = kraw(n, 4, 2, 1, -1)
    nv = n + 1
    rows, rhs = [], []
    # equalities
    r = [fmpq(0)] * nv; r[0] = fmpq(1); rows.append(r); rhs.append(fmpq(1))
    rows.append([fmpq(1)] * nv); rhs.append(fmpq(3 ** n, 2))
    rows.append([fmpq(j) for j in range(nv)]); rhs.append(fmpq(4 * n * 3 ** n, 9))
    for j in range(1, t + 1):
        r = [fmpq(0)] * nv; r[j] = fmpq(1); rows.append(r); rhs.append(fmpq(0))
        rows.append([fmpq(K[j][i]) for i in range(nv)]); rhs.append(fmpq(0))
    neq = len(rows)
    # inequalities with slacks: B_m - A_m - s = 0  ->  (2/3^n) K[m] a - e_m - s = 0
    ineq = []
    for m_ in range(2, n + 1):
        ineq.append([fmpq(2 * K[m_][i], 3 ** n) - (1 if i == m_ else 0) for i in range(nv)])
    for m_ in range(n + 1):
        ineq.append([fmpq(Ks[m_][i]) for i in range(nv)])
    ns = len(ineq)
    A = [row + [fmpq(0)] * ns for row in rows]
    b = list(rhs)
    for k, row in enumerate(ineq):
        A.append(row + [fmpq(-1) if kk == k else fmpq(0) for kk in range(ns)])
        b.append(fmpq(0))
    c = [fmpq((-1) ** (n - j), 2 ** (n - j)) for j in range(nv)] + [fmpq(0)] * ns
    st, val, x = simplex(A, b, c)
    if st != 'optimal':
        return st, None
    return st, fmpq(4, 3 ** n) * val


if __name__ == '__main__':
    n0, n1 = int(sys.argv[1]), int(sys.argv[2])
    for n in range(n0, n1 + 1):
        out = []
        for t in range(1, n // 2 + 1):
            st, q = lp(n, t)
            if st == 'infeasible':
                out.append(f'{t}:inf'); break
            out.append(f'{t}:{float(q):+.4g}')
            if q > 0 and len(out) >= 2 and out[-2].endswith('inf') is False and float(q) > 0:
                pass
        print(f'n={n:2d} ' + ' '.join(out), flush=True)
