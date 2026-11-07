"""Checks for Theorem 7.18 (the n-copy positivity threshold on a real plane is a_n, for every n).

On a real plane P_a compresses to Q_a = (2a - 1/2) Ad_Y o Y_t, where Y_t is the Pauli map with eigenvalues
(1, t, 0, t) on (I, X, Y, Z) and t = 1/(4a - 1).  With u = 2t, a >= a_n is u <= 2^(1/n).  At u = 2^(1/n) the Choi
matrix of Y_t^(x)n is P + P^Gamma with the Bell-diagonal operator P >= 0 given in the theorem.

1. gamma: the action of the partial transpose on Bell-diagonal operators, computed from explicit 4 x 4 matrices.
2. the three one-pair identities (exact, SymPy, u a symbol): gamma f1 = f1, gamma f~1 = u r, gamma r = f~1 / u.
3. the n-copy identity p + gamma^n p = f^n on Bell strings, exact in the number field Q(2^(1/n)), n = 1..4.
4. explicit matrices, n = 1..5 (floating point): the one-copy Choi matrix has Bell eigenvalues (1+u, 1, 1, 1-u)/4,
   P >= 0 and P + P^Gamma = C^(x)n to rounding, with the partial transpose taken directly on the matrix.
5. control: at u = 1.01 * 2^(1/n) the same formula does not give a decomposition (the identity fails).
Ends with PASS or FAIL."""
import itertools, sys
import numpy as np
import sympy as sp

ok = True
# 1. gamma from explicit matrices
bell = [np.array([1, 0, 0, 1]) / np.sqrt(2), np.array([1, 0, 0, -1]) / np.sqrt(2),
        np.array([0, 1, 1, 0]) / np.sqrt(2), np.array([0, 1, -1, 0]) / np.sqrt(2)]
def pt2(M):   # partial transpose of the second qubit of a 4 x 4 matrix
    return M.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)
G = np.array([[np.real(bell[i] @ pt2(np.outer(bell[j], bell[j])) @ bell[i]) for j in range(4)] for i in range(4)])
gamma_expected = 0.5 * np.array([[1, 1, 1, -1], [1, 1, -1, 1], [1, -1, 1, 1], [-1, 1, 1, 1]])
offdiag = max(abs(bell[i] @ pt2(np.outer(bell[j], bell[j])) @ bell[k]) for i in range(4) for j in range(4)
              for k in range(4) if i != k)
g_ok = np.allclose(G, gamma_expected) and offdiag < 1e-12
print(f"1. partial transpose maps Bell projectors to Bell-diagonal operators: {offdiag < 1e-12}; gamma as stated: {g_ok}")
ok &= g_ok
# 2. one-pair identities
u = sp.symbols('u', positive=True)
gam = sp.Rational(1, 2) * sp.Matrix([[1, 1, 1, -1], [1, 1, -1, 1], [1, -1, 1, 1], [-1, 1, 1, 1]])
f1 = sp.Matrix([1 + u, 1, 1, 1 - u]); ft1 = sp.Matrix([1 + u, 1, 1, u - 1]); r = sp.Matrix([2 / u, 1, 1, 0])
ids = [sp.simplify(gam * f1 - f1), sp.simplify(gam * ft1 - u * r), sp.simplify(gam * r - ft1 / u)]
i_ok = all(all(x == 0 for x in v) for v in ids)
print(f"2. gamma f1 = f1, gamma f~1 = u r, gamma r = f~1/u (exact): {i_ok}")
ok &= i_ok
# 3. n-copy identity, exact in the number field Q(w), w = 2^(1/n): elements are coefficient lists of 1, w, ..., w^(n-1)
from fractions import Fraction as Fr
def nf(n):
    def mul(a, b):
        out = [Fr(0)] * n
        for i, x in enumerate(a):
            if x == 0: continue
            for j, y in enumerate(b):
                if y == 0: continue
                k = i + j
                if k >= n: out[k - n] += 2 * x * y        # w^n = 2
                else: out[k] += x * y
        return out
    def add(a, b): return [x + y for x, y in zip(a, b)]
    def sc(c, a): return [c * x for x in a]
    one = [Fr(1)] + [Fr(0)] * (n - 1)
    w = ([Fr(0), Fr(1)] + [Fr(0)] * (n - 2)) if n > 1 else [Fr(2)]
    winv = sc(Fr(1, 2), ([Fr(0)] * (n - 1) + [Fr(1)]) if n > 1 else [Fr(1)])   # w^-1 = w^(n-1)/2
    return mul, add, sc, one, w, winv
for n in range(1, 5):
    mul, add, sc, one, w, winv = nf(n)
    zero = [Fr(0)] * n
    fv = [add(one, w), one, one, add(one, sc(-1, w))]          # 1+u, 1, 1, 1-u
    ftv = [add(one, w), one, one, add(w, sc(-1, one))]         # 1+u, 1, 1, u-1
    rv = [sc(2, winv), one, one, zero]                         # 2/u, 1, 1, 0
    def prodv(vec, J):
        out = one
        for j in J: out = mul(out, vec[j])
        return out
    strings = list(itertools.product(range(4), repeat=n))
    p = {J: add(sc(Fr(1, 2), add(prodv(fv, J), prodv(ftv, J))), sc(-1, prodv(rv, J))) for J in strings}
    g = [[Fr(1, 2) * x for x in row] for row in [[1, 1, 1, -1], [1, 1, -1, 1], [1, -1, 1, 1], [-1, 1, 1, 1]]]
    good = True
    wf = 2 ** (1 / n); minp = None
    for J in strings:
        gp = zero
        for K in strings:
            c = Fr(1)
            for s_ in range(n): c *= g[J[s_]][K[s_]]
            if c != 0: gp = add(gp, sc(c, p[K]))
        res = add(add(p[J], gp), sc(-1, prodv(fv, J)))
        good &= all(x == 0 for x in res)
        val = float(sum(float(x) * wf**i for i, x in enumerate(p[J])))
        minp = val if minp is None else min(minp, val)
    e_ok = good and minp > -1e-12
    print(f"3. n = {n}: p + gamma^n p = f^n exactly in Q(2^(1/n)): {good}; min p = {minp:.6f}")
    ok &= e_ok
# 4. explicit matrices
I2 = np.eye(2); X = np.array([[0, 1], [1, 0]]); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1, -1])
def ptrans_all(M, n):
    T = M.reshape([2] * (4 * n)); perm = list(range(4 * n))
    for s in range(n):
        rb, cb = 2 * s + 1, 2 * n + 2 * s + 1
        perm[rb], perm[cb] = perm[cb], perm[rb]
    return T.transpose(perm).reshape(4**n, 4**n)
for n in range(1, 6):
    uu = 2 ** (1 / n); t = uu / 2
    C1 = 0.25 * sum(l * np.kron(P, P.conj()) for l, P in zip([1, t, 0, t], (I2, X, Y, Z)))
    ev1 = np.array([np.real(b @ C1 @ b) for b in bell]) * 4
    Cn = C1
    for _ in range(n - 1): Cn = np.kron(Cn, C1)
    a1, b1, c1 = np.array([1 + uu, 1, 1, 1 - uu]), np.array([1 + uu, 1, 1, uu - 1]), np.array([2 / uu, 1, 1, 0])
    Pn = np.zeros((4**n, 4**n)); mn = np.inf
    for J in itertools.product(range(4), repeat=n):
        val = 0.5 * (np.prod(a1[list(J)]) + np.prod(b1[list(J)])) - np.prod(c1[list(J)])
        mn = min(mn, val)
        v = bell[J[0]]
        for j in J[1:]: v = np.kron(v, bell[j])
        Pn += val / 4**n * np.outer(v, v)
    err = np.abs(Pn + ptrans_all(Pn, n) - Cn).max()
    m_ok = np.allclose(ev1, [1 + uu, 1, 1, 1 - uu]) and mn > -1e-12 and err < 1e-12
    print(f"4. n = {n}: Bell eigenvalues x4 {np.round(ev1, 6)}; min eigenvalue of 4^n P = {mn:.2e}; "
          f"max |P + P^Gamma - C^n| = {err:.1e}")
    ok &= m_ok
# 5. control above the threshold
n = 3; uu = 1.01 * 2 ** (1 / 3)
a1, b1, c1 = np.array([1 + uu, 1, 1, 1 - uu]), np.array([1 + uu, 1, 1, uu - 1]), np.array([2 / uu, 1, 1, 0])
gm = gamma_expected
res = 0
for J in itertools.product(range(4), repeat=n):
    p = lambda K: 0.5 * (np.prod(a1[list(K)]) + np.prod(b1[list(K)])) - np.prod(c1[list(K)])
    gp = sum(np.prod([gm[J[s], K[s]] for s in range(n)]) * p(K) for K in itertools.product(range(4), repeat=n))
    res = max(res, abs(p(J) + gp - np.prod(a1[list(J)])))
c_ok = res > 1e-3
print(f"5. control at u = 1.01 * 2^(1/3): the identity fails (max residual {res:.3f}): {c_ok}")
ok &= c_ok
print("PASS" if ok else "FAIL")
