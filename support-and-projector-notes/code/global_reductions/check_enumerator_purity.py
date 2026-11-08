"""Section 7, purity two: the enumerator relations do not imply q_n >= 0 for rank-two projectors whose range is
pure two-uniform (A_1 = A_2 = B_1 = B_2 = 0).  A floating-point LP (scipy/HiGHS) proposes an optimal vertex; the vertex is
then recomputed in exact rational arithmetic from its active constraints, and every relation is re-checked exactly.
Normalisation as in Section 7, for a rank-K projector: A_j = K^-2 sum_{wt E=j} |Tr Pi E|^2, B_j = K^-1 sum Tr(Pi E Pi E^*),
A_0 = B_0 = 1, sum_j A_j = 3^n/K, B_j >= A_j >= 0, B(x,y) = (K/3^n) A(x+8y, x-y), coefficients of A(4x+2y, x-y) >= 0,
q_n(Pi) = (K^2/3^n) sum_j (-1/2)^{n-j} A_j."""
from fractions import Fraction as Fr
from math import comb
import numpy as np
import sympy as sp
from scipy.optimize import linprog


def coeffs(n, p, q, r, s):
    """M[m][j] = coefficient of x^{n-m} y^m in (p x + q y)^{n-j} (r x + s y)^j."""
    M = [[0] * (n + 1) for _ in range(n + 1)]
    for j in range(n + 1):
        for h in range(n - j + 1):
            for g in range(j + 1):
                M[h + g][j] += comb(n - j, h) * p ** (n - j - h) * q ** h * comb(j, g) * r ** (j - g) * s ** g
    return M


def relations(n, K, t):
    """list of (row, kind, rhs) on the variables A_0..A_n, all exact."""
    Bm = coeffs(n, 1, 8, 1, -1)          # B_m = (K/3^n) sum_j Bm[m][j] A_j
    Sm = coeffs(n, 4, 2, 1, -1)
    e = lambda j: [Fr(int(i == j)) for i in range(n + 1)]
    rows = [(e(0), 'eq', Fr(1)), ([Fr(1)] * (n + 1), 'eq', Fr(3 ** n, K))]
    for m in range(1, n + 1):
        Brow = [Fr(K, 3 ** n) * c for c in Bm[m]]
        if m <= t:
            rows += [(e(m), 'eq', Fr(0)), (Brow, 'eq', Fr(0))]
        else:
            rows += [([b - a for b, a in zip(Brow, e(m))], 'ge', Fr(0))]      # B_m >= A_m
    rows += [([Fr(c) for c in Sm[m]], 'ge', Fr(0)) for m in range(n + 1)]
    rows += [(e(j), 'ge', Fr(0)) for j in range(n + 1)]
    return rows


def exact_vertex(n, K, t):
    rows = relations(n, K, t)
    obj = [Fr(K * K, 3 ** n) * Fr(-1, 2) ** (n - j) for j in range(n + 1)]
    sc = lambda r: max(abs(float(c)) for c in r) or 1.0
    eq = [(r, rhs) for r, k, rhs in rows if k == 'eq']
    ge = [r for r, k, rhs in rows if k == 'ge']
    res = linprog([float(c) for c in obj],
                  A_ub=np.array([[-float(c) / sc(r) for c in r] for r in ge]), b_ub=np.zeros(len(ge)),
                  A_eq=np.array([[float(c) / sc(r) for c in r] for r, _ in eq]),
                  b_eq=np.array([float(rhs) / sc(r) for r, rhs in eq]), bounds=[(None, None)] * (n + 1), method='highs')
    x = res.x
    slack = lambda r, rhs: abs(sum(float(c) * xi for c, xi in zip(r, x)) - float(rhs)) / sc(r)
    order = sorted(rows, key=lambda z: 0.0 if z[1] == 'eq' else slack(z[0], z[2]))
    Msel, rhs = [], []
    for r, k, h in order:
        if sp.Matrix(Msel + [r]).rank() > len(Msel):
            Msel.append(r); rhs.append(h)
        if len(Msel) == n + 1:
            break
    sol = sp.Matrix(Msel).LUsolve(sp.Matrix(rhs))
    A = [Fr(int(sp.fraction(v)[0]), int(sp.fraction(v)[1])) for v in sol]
    feasible = all((sum(c * a for c, a in zip(r, A)) == h) if k == 'eq' else (sum(c * a for c, a in zip(r, A)) >= h)
                   for r, k, h in rows)
    qv = sum(c * a for c, a in zip(obj, A))
    return feasible, qv, A, res.fun


ok = True
for n, K, t in ((10, 2, 2), (11, 2, 2), (10, 3, 2), (11, 3, 2)):
    feas, qv, A, fl = exact_vertex(n, K, t)
    print(f"n = {n}, rank {K}, purity {t}: exact feasible enumerator {feas}; q_n = {qv} = {float(qv):.6f}; "
          f"q_n/||Pi||^2 = {float(qv) / K:.4f} (LP optimum, floating point: {fl:.6f})")
    if K == 2 and n == 10:
        print('   A =', [str(a) for a in A])
    ok &= feas and qv < 0 and abs(float(qv) - fl) < 1e-9
print('purity-two enumerator relations allow q_n < 0 (exact):', ok)
