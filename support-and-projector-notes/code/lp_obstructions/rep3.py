"""Three-replica twirl algebra for rank-two projectors on n qutrit copies.

Z = Pi (x) Pi (x) Pi on three replicas; T = twirl of Z over prod_nu U_nu^{(x)3} lies in (x)_nu C[S_3].
C[S_3] acting on (C^3)^{(x)3} = (+)_lam V_lam (x) S_lam, lam in {t (trivial, dim V 10), s (standard, dim V 8, dim S 2),
a (sign, dim V 1)}. After averaging over copy permutations, the block of T with (a, m, b) copies in (t, s, a) is
(+)_j X_j (x) I_{M_j} on (C^2)^{(x)m}, X_j acting on the spin-j space V_j = Sym^{2j} (x) det^{m/2-j}.
"""
import itertools
from math import comb, factorial
import numpy as np

PERMS = list(itertools.permutations(range(3)))          # p maps i -> p[i]
E = (0, 1, 2)


def compose(p, q):
    return tuple(p[q[i]] for i in range(3))


def inverse(p):
    r = [0] * 3
    for i, pi in enumerate(p):
        r[pi] = i
    return tuple(r)


def ncycles(p):
    seen, c = set(), 0
    for i in range(3):
        if i not in seen:
            c += 1
            j = i
            while j not in seen:
                seen.add(j); j = p[j]
    return c


def sign(p):
    return (-1) ** (3 - ncycles(p))


def permmat(p):
    M = np.zeros((3, 3))
    for k in range(3):
        M[p[k], k] = 1.0
    return M


_F = np.array([[1, 1], [-1, 1], [0, -2]], float) / np.array([np.sqrt(2), np.sqrt(6)])
STD = {p: _F.T @ permmat(p) @ _F for p in PERMS}           # orthogonal 2-dim irrep, homomorphism
for p in PERMS:
    for q in PERMS:
        assert np.allclose(STD[p] @ STD[q], STD[compose(p, q)])

DIMV = {'t': 10, 's': 8, 'a': 1}


def dimM(m, j2):
    """dimension of the S_m irrep paired with spin j (j2 = 2j) in (C^2)^{(x)m}."""
    k = (m - j2) // 2
    return comb(m, k) - (comb(m, k - 1) if k >= 1 else 0)


# ---- polynomials in y as numpy coefficient arrays (index = power of y) ----
def padd(a, b):
    n = max(len(a), len(b)); r = np.zeros(n)
    r[:len(a)] += a; r[:len(b)] += b
    return r


def pmul(a, b):
    return np.convolve(a, b)


def ppow(a, k):
    r = np.array([1.0])
    for _ in range(k):
        r = pmul(r, a)
    return r


def sym_power_poly(A, d):
    """Sym^d of a 2x2 matrix of polynomials, in the orthonormal basis of symmetric tensors.
    A[r][c] are polynomial arrays. Returns (d+1)x(d+1) list of polynomial arrays."""
    # monomial basis e_i = u^{d-i} v^i ; A maps u -> A00 u + A10 v, v -> A01 u + A11 v
    S = [[np.array([0.0]) for _ in range(d + 1)] for _ in range(d + 1)]
    for i in range(d + 1):
        # (A00 u + A10 v)^{d-i} (A01 u + A11 v)^i
        for s in range(d - i + 1):           # choose v from first factor s times
            c1 = comb(d - i, s)
            t1 = pmul(ppow(A[0][0], d - i - s), ppow(A[1][0], s))
            for r in range(i + 1):           # choose v from second factor r times
                c2 = comb(i, r)
                t2 = pmul(ppow(A[0][1], i - r), ppow(A[1][1], r))
                row = s + r                  # power of v
                S[row][i] = padd(S[row][i], c1 * c2 * pmul(t1, t2))
    # orthonormal: S_on[r][c] = S[r][c] * sqrt(C(d,r)/C(d,c))
    for r in range(d + 1):
        for c in range(d + 1):
            S[r][c] = S[r][c] * np.sqrt(comb(d, r) / comb(d, c))
    return S


def pi_j(A, m, j2):
    """action of A^{(x)m} on the spin-j isotypic space: Sym^{2j}(A) det(A)^{m/2 - j}."""
    S = sym_power_poly(A, j2)
    det = padd(pmul(A[0][0], A[1][1]), -pmul(A[0][1], A[1][0]))
    dp = ppow(det, (m - j2) // 2)
    return [[pmul(S[r][c], dp) for c in range(j2 + 1)] for r in range(j2 + 1)]


def blocks(n):
    out = []
    for m in range(n + 1):
        for a in range(n - m + 1):
            b = n - m - a
            for j2 in range(m % 2, m + 1, 2):
                out.append((a, m, b, j2))
    return out


def element_images(f):
    """f: dict perm -> polynomial array. Returns (f_t, f_a, f_s) with f_s a 2x2 polynomial matrix."""
    ft = np.array([0.0]); fa = np.array([0.0])
    fs = [[np.array([0.0]), np.array([0.0])], [np.array([0.0]), np.array([0.0])]]
    for p, c in f.items():
        ft = padd(ft, c); fa = padd(fa, sign(p) * c)
        for r in range(2):
            for s in range(2):
                fs[r][s] = padd(fs[r][s], STD[p][r, s] * c)
    return ft, fa, fs


def product_functional(n, f):
    """Coefficient matrices: for each block (a,m,b,j2) a list over powers of y of (j2+1)x(j2+1) matrices M_k such
    that Tr[T f^{(x)n}] = sum_blocks sum_k y^k Tr[X M_k]."""
    ft, fa, fs = element_images(f)
    res = {}
    for (a, m, b, j2) in blocks(n):
        mult = factorial(n) // (factorial(a) * factorial(m) * factorial(b)) * DIMV['t'] ** a * DIMV['s'] ** m
        scal = mult * dimM(m, j2) * pmul(ppow(ft, a), ppow(fa, b))
        P = pi_j(fs, m, j2)
        deg = max(len(pmul(scal, P[r][c])) for r in range(j2 + 1) for c in range(j2 + 1))
        mats = np.zeros((deg, j2 + 1, j2 + 1))
        for r in range(j2 + 1):
            for c in range(j2 + 1):
                poly = pmul(scal, P[r][c])
                mats[:len(poly), r, c] = poly
        # Tr[X P] = sum_{r,c} X[c,r] P[r,c]; store transposed so that functional = sum X * M
        res[(a, m, b, j2)] = np.transpose(mats, (0, 2, 1))
    return res


def const(c):
    return np.array([float(c)])


def poly(c0, c1):
    return np.array([float(c0), float(c1)])


T12 = (1, 0, 2); T13 = (2, 1, 0); T23 = (0, 2, 1); C1 = (1, 2, 0); C2 = (2, 0, 1)


def sign_projector(m, j2, b):
    """(1/6) sum_tau sgn(tau)^{1+b} pi_j(rho_std(tau)) on the spin-j space."""
    P = np.zeros((j2 + 1, j2 + 1))
    for p in PERMS:
        A = [[const(STD[p][r, s]) for s in range(2)] for r in range(2)]
        M = pi_j(A, m, j2)
        P += sign(p) ** (1 + b) * np.array([[M[r][c][0] for c in range(j2 + 1)] for r in range(j2 + 1)]) / 6
    return P


def group_action(m, j2, p):
    A = [[const(STD[p][r, s]) for s in range(2)] for r in range(2)]
    M = pi_j(A, m, j2)
    return np.array([[M[r][c][0] for c in range(j2 + 1)] for r in range(j2 + 1)])
