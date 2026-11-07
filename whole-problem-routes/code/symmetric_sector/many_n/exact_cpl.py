"""Exact rational verification of word-span certificates in the coupled basis.

Every matrix element of a.b and b^+.a in the orthonormal coupled basis (cpl.py) is sign * sqrt(r) with r rational
(6j^2 times integer reduced matrix elements squared).  Rescale each basis state i by sqrt(s_i), s_i a square-free
integer chosen so that r * s_i * s_i' is a rational square for every nonzero entry (square classes propagated along the
nonzero entries; consistency is asserted).  In the rescaled bases all operators are rational matrices:
    O' = S'^{-1/2} O S^{1/2},
and the quadratic form of Q = sum_t Xt_t^T G_t Xt_t on the (n,n) spin-L block in the rescaled basis is
    Q'' = sum_t Xt'_t^T S_m G'_t Xt'_t   (congruent to the orthonormal block by S_n^{1/2}),
with G_t = g_t(C_m) a polynomial in C_m = E_ab E_ba (eigenvalues (m-j)(m-j+1)) interpolating g_t(j) = sum_s phi omega_j(s).
Q'' is then tested by exact LDL^T over Q (python-flint fmpq).  Nothing here uses floating point.
"""
import math
from fractions import Fraction
from functools import lru_cache
import flint
from cpl import basis, sixj_exact, omega


def sqfree(x):
    """square-free part of a positive integer (trial division; inputs here are products of small factors)."""
    x = int(x)
    assert x > 0
    out = 1
    p = 2
    while p * p <= x:
        e = 0
        while x % p == 0:
            x //= p; e += 1
        if e % 2:
            out *= p
        p += 1 if p == 2 else 2
    return out * x


def sq_class(fr):
    fr = Fraction(fr)
    return sqfree(fr.numerator * fr.denominator)


def rat_sqrt(fr):
    fr = Fraction(fr)
    a, b = fr.numerator, fr.denominator
    ra, rb = math.isqrt(a), math.isqrt(b)
    assert ra * ra == a and rb * rb == b, ('not a square', fr)
    return Fraction(ra, rb)


def red2_cr(N, l, lp):
    if lp == l + 1:
        return (1, (l + 1) * (N + l + 3))
    if lp == l - 1:
        return (-1, l * (N - l + 2))
    return (0, 0)


def red2_an(N, l, lp):
    if lp == l + 1:
        return (1, (l + 1) * (N - l)) if N - l > 0 else (0, 0)
    if lp == l - 1:
        return (-1, l * (N + l + 1))
    return (0, 0)


def entries_sq(n1, n2, L, k1, k2):
    """Nonzero entries of T.U from sector (n1,n2) to (n1+-1, n2+-1), spin L, as {(i', i): (sign, r)}, value = sign*sqrt(r)."""
    n1p = n1 + (1 if k1 == 'c' else -1); n2p = n2 + (1 if k2 == 'c' else -1)
    B = basis(n1, n2, L); Bp = basis(n1p, n2p, L)
    idx = {b: i for i, b in enumerate(Bp)}
    f1 = red2_cr if k1 == 'c' else red2_an
    f2 = red2_cr if k2 == 'c' else red2_an
    out = {}
    for col, (j1, j2) in enumerate(B):
        for j1p in (j1 - 1, j1 + 1):
            for j2p in (j2 - 1, j2 + 1):
                i = idx.get((j1p, j2p))
                if i is None:
                    continue
                s6, q6 = sixj_exact(j1p, j2p, L, j2, j1, 1)
                s1, q1 = f1(n1, j1, j1p); s2, q2 = f2(n2, j2, j2p)
                if s6 == 0 or s1 == 0 or s2 == 0:
                    continue
                ph = -1 if (j1 + j2p + L) % 2 else 1
                out[(i, col)] = (ph * s6 * s1 * s2, q6 * q1 * q2)
    return out


class ExactL:
    """All rescaled rational data for one spin L and top level n."""

    def __init__(self, n, L):
        self.n, self.L = n, L
        # sectors: ('d', m) = (m, m); ('o', m) = (m-1, m+1)
        self.sec_basis = {}
        for m in range(0, n + 1):
            self.sec_basis[('d', m)] = basis(m, m, L)
            if m >= 1:
                self.sec_basis[('o', m)] = basis(m - 1, m + 1, L)
        # edges
        self.ab_e = {m: entries_sq(m, m, L, 'a', 'a') for m in range(1, n + 1)}      # ('d',m) -> ('d',m-1)
        self.eba_e = {m: entries_sq(m, m, L, 'a', 'c') for m in range(1, n + 1)}     # ('d',m) -> ('o',m)
        self._assign_classes()

    def _assign_classes(self):
        s = {key: [None] * len(b) for key, b in self.sec_basis.items()}
        adj = {}
        def add_edge(k1, i1, k2, i2, r):
            c = sq_class(r)
            adj.setdefault((k1, i1), []).append((k2, i2, c))
            adj.setdefault((k2, i2), []).append((k1, i1, c))
        for m, E in self.ab_e.items():
            for (ip, i), (sg, r) in E.items():
                add_edge(('d', m), i, ('d', m - 1), ip, r)
        for m, E in self.eba_e.items():
            for (ip, i), (sg, r) in E.items():
                add_edge(('d', m), i, ('o', m), ip, r)
        for key, b in self.sec_basis.items():
            for i in range(len(b)):
                if s[key][i] is not None:
                    continue
                s[key][i] = 1
                stack = [(key, i)]
                while stack:
                    k0, i0 = stack.pop()
                    for (k1, i1, c) in adj.get((k0, i0), []):
                        want = sqfree(c * s[k0][i0])
                        if s[k1][i1] is None:
                            s[k1][i1] = want; stack.append((k1, i1))
                        else:
                            assert s[k1][i1] == want, 'square classes inconsistent'
        self.s = s

    def rescaled(self, E, src, dst):
        """Rational matrix O' = S_dst^{-1/2} O S_src^{1/2} as fmpq_mat."""
        nr = len(self.sec_basis[dst]); nc = len(self.sec_basis[src])
        M = [[0] * nc for _ in range(nr)]
        for (ip, i), (sg, r) in E.items():
            M[ip][i] = sg * rat_sqrt(Fraction(r) * self.s[src][i] / self.s[dst][ip])
        return flint.fmpq_mat(nr, nc, [flint.fmpq(x.numerator, x.denominator) if isinstance(x, Fraction) else flint.fmpq(x)
                                        for row in M for x in row])

    def Sdiag(self, key):
        v = self.s[key]
        d = len(v)
        return flint.fmpq_mat(d, d, [flint.fmpq(v[i]) if i == j else flint.fmpq(0) for i in range(d) for j in range(d)])

    def Sinv(self, key):
        v = self.s[key]
        d = len(v)
        return flint.fmpq_mat(d, d, [flint.fmpq(1, v[i]) if i == j else flint.fmpq(0) for i in range(d) for j in range(d)])


def interp_coeffs(xs, ys):
    """Exact monomial coefficients of the interpolating polynomial through (xs, ys) (Newton form expanded)."""
    nn = len(xs)
    dd = [Fraction(y) for y in ys]
    coef = [dd[0]]
    for j in range(1, nn):
        for i in range(nn - 1, j - 1, -1):
            dd[i] = (dd[i] - dd[i - 1]) / (xs[i] - xs[i - j])
        coef.append(dd[j])
    # expand Newton form: p(x) = coef0 + coef1 (x - x0) + coef2 (x - x0)(x - x1) + ...
    poly = [Fraction(0)] * nn
    basisp = [Fraction(1)]
    for j in range(nn):
        for k, b in enumerate(basisp):
            poly[k] += coef[j] * b
        if j < nn - 1:
            nb = [Fraction(0)] * (len(basisp) + 1)
            for k, b in enumerate(basisp):
                nb[k + 1] += b
                nb[k] -= xs[j] * b
            basisp = nb
    return poly


def q_block(EL, phi):
    """phi: {(k,s,t): Fraction} coefficients of the monomial words M(k,s,t) (symbol p^k q^s r^t).
    Returns the exact rational matrix Q'' (rescaled basis of the (n,n) spin-L block)."""
    n, L = EL.n, EL.L
    dn = len(EL.sec_basis[('d', n)])
    Q = flint.fmpq_mat(dn, dn)
    # Xt'_t: rescaled (m!/n!) (a.b)^t
    X = flint.fmpq_mat(dn, dn, [flint.fmpq(1) if i == j else flint.fmpq(0) for i in range(dn) for j in range(dn)])
    bydt = {}
    for (k, s, t), v in phi.items():
        bydt.setdefault(t, {})[s] = v
    for t in range(0, n + 1):
        m = n - t
        if t > 0:
            A = EL.rescaled(EL.ab_e[m + 1], ('d', m + 1), ('d', m))
            X = (A * X) * flint.fmpq(1, m + 1)
        dm = len(EL.sec_basis[('d', m)])
        if dm == 0 or t not in bydt:
            continue
        # g_t(j) = sum_s phi(m-s, s, t) omega_j^{(m)}(s), interpolated in C_m with eigenvalues (m-j)(m-j+1)
        xs = [Fraction((m - j) * (m - j + 1)) for j in range(m + 1)]
        ys = [sum(Fraction(v) * omega(m, j, s) for s, v in bydt[t].items()) for j in range(m + 1)]
        if m >= 1:
            E = EL.rescaled(EL.eba_e[m], ('d', m), ('o', m))
            C = EL.Sinv(('d', m)) * (E.transpose() * (EL.Sdiag(('o', m)) * E))
        else:
            C = flint.fmpq_mat(dm, dm)
        coef = interp_coeffs(xs, ys)
        G = flint.fmpq_mat(dm, dm)
        Id = flint.fmpq_mat(dm, dm, [flint.fmpq(1) if i == j else flint.fmpq(0) for i in range(dm) for j in range(dm)])
        for c in reversed(coef):
            G = C * G + Id * flint.fmpq(c.numerator, c.denominator)
        Q += X.transpose() * ((EL.Sdiag(('d', m)) * G) * X)
    return Q


def ldl_psd(M):
    """Exact LDL^T PSD test of a symmetric fmpq_mat. Returns (is_psd, rank, n_zero)."""
    n = M.nrows()
    A = [[M[i, j] for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(i):
            assert A[i][j] == A[j][i], 'not symmetric'
    rank = 0; zeros = 0
    act = list(range(n))
    while act:
        # pivot: largest diagonal
        pv = max(act, key=lambda i: A[i][i])
        d = A[pv][pv]
        if d < 0:
            return False, rank, zeros
        if d == 0:
            for i in act:
                for j in act:
                    if A[i][j] != 0:
                        return False, rank, zeros
            zeros += len(act)
            break
        act.remove(pv); rank += 1
        rp = A[pv]
        for i in act:
            f = A[i][pv] / d
            if f == 0:
                continue
            Ai = A[i]
            for j in act:
                if rp[j] != 0:
                    Ai[j] -= f * rp[j]
    return True, rank, zeros
