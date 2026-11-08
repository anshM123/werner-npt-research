"""Checks for Proposition M' (balanced scalar form of the map route).  P_a(X) = a Tr(X) I_3 - (X + X^T)/4, m = a - 1/2.
(1) exact (sympy, polynomials in a): the flag kernel g(p,q,r,s) = h_r^T P_a^{(x)2}(h_p h_q^T) h_s of the flags
    h_0 = (|00>+|11>+|22>)/sqrt3, h_1 = (|01>+|12>+|20>)/sqrt3 is A = a^2 - a/3 + 1/6 (all equal), B = a^2 - a/3
    (p = q, r = s, p != r), K = (5 - 8a)/48 (the other two-plus-two patterns) and 0 (odd patterns);
(2) exact: s_r(g_r) = [m^r + 2a^r + 4(-1/4)^r]/3 for g_r = (|0^r> + |1^r> + |2^r>)/sqrt3, r = 2..5;
(3) floating point: the two-flag identity s_{n+k}(z_k) = m^k[s_n(x)+s_n(y)] + 2a^k b_n(x,y) + 4(-1/4)^k kappa_n(x,y);
(4) exact (Fractions): at a = 7/10 the plane GHZ vector x = (p^{(x)7} + conj(p)^{(x)7})/sqrt2, p = (e_0 + i e_1)/sqrt2, has
    s_7(x) = -217031/1280000000, through the compression Q_a(Y) = a Tr(Y) I_2 - (Y + Y^T)/4 to span{e_0, e_1};
(5) exact: E_4(7/10) (8/5)^7 < 217031/1280000000 with E_l = 11 r^{8l} + 22 h^{8l} + 78 r^{4l} h^{8l}, r = B/A, h = |K|/A;
    (7/2)(15/23)^4 < 1 and 444 (7/2)^3 (15/23)^24 < 1; and at a = 3/4, s_3(g_3)/a^3 = 17/27, A/a^2 = 23/27.
usage: python balanced_scalar_check.py"""
import itertools
from fractions import Fraction as Fr
import numpy as np
import sympy as sp

ok = True
a = sp.symbols('a')


def Pmap(X, par):                 # one qutrit, sympy or numpy
    if isinstance(X, np.ndarray):
        return par * np.trace(X) * np.eye(3) - (X + X.T) / 4
    return par * X.trace() * sp.eye(3) - (X + X.T) / 4


def P2(Y, par):                   # P_a (x) P_a on a 9x9 operator (two qutrits), sympy
    T = {}
    out = sp.zeros(9, 9)
    # apply site maps via the basis expansion Y = sum Y[(i1 i2),(j1 j2)] |i1><j1| (x) |i2><j2|
    for i1, i2, j1, j2 in itertools.product(range(3), repeat=4):
        c = Y[3 * i1 + i2, 3 * j1 + j2]
        if c == 0:
            continue
        E1 = sp.zeros(3, 3); E1[i1, j1] = 1
        E2 = sp.zeros(3, 3); E2[i2, j2] = 1
        out += c * sp.kronecker_product(Pmap(E1, par), Pmap(E2, par))
    return out


def svec(entries):
    v = sp.zeros(9, 1)
    for (i, j) in entries:
        v[3 * i + j] = 1
    return v / sp.sqrt(3)


h = [svec([(k, k) for k in range(3)]), svec([(k, (k + 1) % 3) for k in range(3)])]
A, B, K = a ** 2 - a / 3 + sp.Rational(1, 6), a ** 2 - a / 3, (5 - 8 * a) / 48
good = True
for p, q, r, s in itertools.product(range(2), repeat=4):
    val = sp.expand((h[r].T * P2(h[p] * h[q].T, a) * h[s])[0, 0])
    want = A if p == q == r == s else B if (p == q and r == s) else K if (p + q + r + s) % 2 == 0 else 0
    good &= sp.expand(val - want) == 0
ok &= good
print(f"(1) flag kernel equals A, B, K, 0 by pattern, as polynomials in a: {good}")

m = a - sp.Rational(1, 2)
good = True
for rr in range(2, 6):
    D = 3 ** rr
    gvec = np.zeros(D)
    for k in range(3):
        gvec[int(str(k) * rr, 3)] = 1 / np.sqrt(3)
    for av in (0.6, 0.75, 0.9):
        X = np.outer(gvec, gvec).reshape((3,) * (2 * rr))
        for site in range(rr):           # apply P_a at each site
            X = np.moveaxis(X, [site, rr + site], [0, 1])
            Xs = X.reshape(3, 3, -1)
            tr = np.trace(Xs, axis1=0, axis2=1)
            Xs = av * np.einsum('ij,k->ijk', np.eye(3), tr) - (Xs + Xs.transpose(1, 0, 2)) / 4
            X = np.moveaxis(Xs.reshape(X.shape), [0, 1], [site, rr + site])
        val = gvec @ X.reshape(D, D) @ gvec
        want = float(((m ** rr + 2 * a ** rr + 4 * sp.Rational(-1, 4) ** rr) / 3).subs(a, av))
        good &= abs(val - want) < 1e-12
ok &= good
print(f"(2) s_r(g_r) = [m^r + 2a^r + 4(-1/4)^r]/3 for r = 2..5 at a = 0.6, 0.75, 0.9: {good}")


def apply_P(X, n, av):
    X = X.reshape((3,) * (2 * n))
    for site in range(n):
        X = np.moveaxis(X, [site, n + site], [0, 1])
        Xs = X.reshape(3, 3, -1)
        tr = np.trace(Xs, axis1=0, axis2=1)
        Xs = av * np.einsum('ij,k->ijk', np.eye(3), tr) - (Xs + Xs.transpose(1, 0, 2)) / 4
        X = np.moveaxis(Xs.reshape(X.shape), [0, 1], [site, n + site])
    return X.reshape(3 ** n, 3 ** n)


rng = np.random.default_rng(5)
worst = 0.0
for n, k in [(1, 1), (1, 2), (2, 1)]:
    for av in (0.6, 0.8):
        mv = av - 0.5
        x = rng.normal(size=3 ** n); y = rng.normal(size=3 ** n)
        f0 = np.zeros(3 ** k); f0[0] = 1
        f1 = np.zeros(3 ** k); f1[int('1' * k, 3)] = 1
        z = np.kron(x, f0) + np.kron(y, f1)
        s = lambda v, N: v @ apply_P(np.outer(v, v), N, av) @ v
        b = y @ apply_P(np.outer(x, x), n, av) @ y
        kap = x @ apply_P(np.outer(x, y), n, av) @ y
        lhs = s(z, n + k)
        rhs = mv ** k * (s(x, n) + s(y, n)) + 2 * av ** k * b + 4 * (-0.25) ** k * kap
        worst = max(worst, abs(lhs - rhs))
ok &= worst < 1e-10
print(f"(3) two-flag identity: max error {worst:.1e}")

# (4) exact plane-GHZ value at a = 7/10 through the compression Q_a on (R^2)^(x)7
av = Fr(7, 10)
n = 7
x = {}
for w in itertools.product(range(2), repeat=n):
    if sum(w) % 2 == 0:
        x[w] = Fr((-1) ** (sum(w) // 2), 8)
# operator X = x x^T as dict over (row word, col word); apply Q_a site by site
X = {(u, v): x[u] * x[v] for u in x for v in x}
for site in range(n):
    Y = {}
    traces = {}
    for (u, v), c in X.items():
        if u[site] == v[site]:
            key = (u[:site] + u[site + 1:], v[:site] + v[site + 1:])
            traces[key] = traces.get(key, 0) + c
    for (ur, vr), c in traces.items():          # a Tr(.) I_2 part
        for b in range(2):
            u = ur[:site] + (b,) + ur[site:]; v = vr[:site] + (b,) + vr[site:]
            Y[(u, v)] = Y.get((u, v), 0) + av * c
    for (u, v), c in X.items():                 # -(Y + Y^T)/4 part at this site
        Y[(u, v)] = Y.get((u, v), 0) - c / 4
        u2 = u[:site] + (v[site],) + u[site + 1:]; v2 = v[:site] + (u[site],) + v[site + 1:]
        Y[(u2, v2)] = Y.get((u2, v2), 0) - c / 4
    X = Y
s7 = sum(x[u] * c * x[v] for (u, v), c in X.items() if u in x and v in x)
good = s7 == Fr(-217031, 1280000000)
ok &= good
print(f"(4) a = 7/10, plane GHZ on 7 copies: s_7(x) = {s7} (expected -217031/1280000000): {good}")

# (5) exact inequalities
Av, Bv, Kv = av ** 2 - av / 3 + Fr(1, 6), av ** 2 - av / 3, (5 - 8 * av) / 48
r_, h_ = Bv / Av, abs(Kv) / Av
E4 = 11 * r_ ** 32 + 22 * h_ ** 32 + 78 * r_ ** 16 * h_ ** 32
c1 = E4 * Fr(8, 5) ** 7 < Fr(217031, 1280000000)
c2 = Fr(7, 2) * Fr(15, 23) ** 4 < 1
c3 = 444 * Fr(7, 2) ** 3 * Fr(15, 23) ** 24 < 1
a34 = Fr(3, 4); m34 = a34 - Fr(1, 2)
s3g3 = (m34 ** 3 + 2 * a34 ** 3 + 4 * Fr(-1, 4) ** 3) / 3
c4 = s3g3 / a34 ** 3 == Fr(17, 27) and (a34 ** 2 - a34 / 3 + Fr(1, 6)) / a34 ** 2 == Fr(23, 27)
c5 = (Av, Bv, Kv) == (Fr(127, 300), Fr(77, 300), Fr(-1, 80))
ok &= c1 and c2 and c3 and c4 and c5
print(f"(5) A, B, K at 7/10 = 127/300, 77/300, -1/80: {c5}; E_4(7/10)(8/5)^7 = {float(E4 * Fr(8, 5) ** 7):.6e} < 1.69555e-4: {c1}; "
      f"(7/2)(15/23)^4 < 1: {c2}; 444(7/2)^3(15/23)^24 < 1: {c3}; constants at 3/4 (17/27, 23/27): {c4}")
print("PASS" if ok else "FAIL")
