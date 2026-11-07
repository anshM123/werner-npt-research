"""Coupled (Racah) basis for O(3)-invariant operators on Sym^n(R^3) (x) Sym^n(R^3), Fock picture.

Modes a_i (inputs) and b_i (outputs), i = 1, 2, 3.  The (n1, n2)-particle sector is Sym^{n1} (x) Sym^{n2}.
Sym^N(R^3) = (+)_{l = N, N-2, ...} H_l (harmonic degree l, multiplicity one), with orthonormal states |N l m> (3D oscillator).
Spin-L multiplicity space of the (n1, n2) sector: orthonormal basis |n1 j1, n2 j2; L> labelled by (j1, j2),
j1 = n1 mod 2, j2 = n2 mod 2, |j1 - j2| <= L <= j1 + j2.

Reduced matrix elements (Edmonds convention <j'm'|T_q|jm> = <j m 1 q|j' m'> <j'||T||j> / sqrt(2j'+1)):
  <N+1, l+1 || a^+ || N, l> = sqrt((l+1)(N+l+3)),  <N+1, l-1 || a^+ || N, l> = -sqrt(l(N-l+2)),
  <N-1, l+1 || a   || N, l> = sqrt((l+1)(N-l)),    <N-1, l-1 || a   || N, l> = -sqrt(l(N+l+1)).
Scalar products of commuting vector operators (Edmonds 7.1.6):
  <j1' j2' L| T.U |j1 j2 L> = (-1)^{j1 + j2' + L} {j1' j2' L; j2 j1 1} <j1'||T||j1> <j2'||U||j2>.
Everything here is checked numerically against brute force in check_cpl.py.
"""
import math
from fractions import Fraction
from functools import lru_cache
import numpy as np


@lru_cache(maxsize=None)
def fact(k):
    return math.factorial(k)


@lru_cache(maxsize=None)
def sixj_exact(a, b, c, d, e, f):
    """6j symbol {a b c; d e f}, integer arguments. Returns (sign, square as Fraction)."""
    def tri(x, y, z):
        return x + y >= z and y + z >= x and z + x >= y
    if min(a, b, c, d, e, f) < 0 or not (tri(a, b, c) and tri(a, e, f) and tri(d, b, f) and tri(d, e, c)):
        return (0, Fraction(0))

    def d2(x, y, z):
        return Fraction(fact(x + y - z) * fact(x - y + z) * fact(-x + y + z), fact(x + y + z + 1))
    pref = d2(a, b, c) * d2(a, e, f) * d2(d, b, f) * d2(d, e, c)
    tmin = max(a + b + c, a + e + f, d + b + f, d + e + c)
    tmax = min(a + b + d + e, b + c + e + f, c + a + f + d)
    s = 0
    for t in range(tmin, tmax + 1):
        s += Fraction((-1) ** t * fact(t + 1),
                      fact(t - a - b - c) * fact(t - a - e - f) * fact(t - d - b - f) * fact(t - d - e - c)
                      * fact(a + b + d + e - t) * fact(b + c + e + f - t) * fact(c + a + f + d - t))
    if s == 0:
        return (0, Fraction(0))
    return (1 if s > 0 else -1, pref * s * s)


def sixj(a, b, c, d, e, f):
    sg, sq = sixj_exact(a, b, c, d, e, f)
    return 0.0 if sg == 0 else sg * math.sqrt(sq.numerator / sq.denominator)


def red_cr(N, l, lp):
    """<N+1, lp || a^+ || N, l>"""
    if lp == l + 1:
        return math.sqrt((l + 1) * (N + l + 3))
    if lp == l - 1:
        return -math.sqrt(l * (N - l + 2))
    return 0.0


def red_an(N, l, lp):
    """<N-1, lp || a || N, l>"""
    if lp == l + 1:
        return math.sqrt((l + 1) * (N - l)) if N - l >= 0 else 0.0
    if lp == l - 1:
        return -math.sqrt(l * (N + l + 1))
    return 0.0


@lru_cache(maxsize=None)
def basis(n1, n2, L):
    out = []
    if n1 < 0 or n2 < 0:
        return tuple(out)
    for j1 in range(n1 % 2, n1 + 1, 2):
        for j2 in range(n2 % 2, n2 + 1, 2):
            if abs(j1 - j2) <= L <= j1 + j2:
                out.append((j1, j2))
    return tuple(out)


def mult(n, L):
    return len(basis(n, n, L))


def scal(n1, n2, L, k1, k2):
    """Matrix of T.U (T = a^+ if k1 == 'c' else a, on the inputs; U likewise on the outputs) from the spin-L space of
    sector (n1, n2) to that of (n1 +- 1, n2 +- 1)."""
    n1p = n1 + (1 if k1 == 'c' else -1)
    n2p = n2 + (1 if k2 == 'c' else -1)
    B = basis(n1, n2, L)
    Bp = basis(n1p, n2p, L)
    idx = {b: i for i, b in enumerate(Bp)}
    M = np.zeros((len(Bp), len(B)))
    f1 = red_cr if k1 == 'c' else red_an
    f2 = red_cr if k2 == 'c' else red_an
    for col, (j1, j2) in enumerate(B):
        for j1p in (j1 - 1, j1 + 1):
            for j2p in (j2 - 1, j2 + 1):
                i = idx.get((j1p, j2p))
                if i is None:
                    continue
                ph = -1.0 if (j1 + j2p + L) % 2 else 1.0
                M[i, col] += ph * sixj(j1p, j2p, L, j2, j1, 1) * f1(n1, j1, j1p) * f2(n2, j2, j2p)
    return M


def ab(n, L):
    """a.b = sum_i a_i b_i : (n, n) -> (n-1, n-1), spin L."""
    return scal(n, n, L, 'a', 'a')


def Eba(m, L):
    """b^+.a : (m, m) -> (m-1, m+1), spin L."""
    return scal(m, m, L, 'a', 'c')


def swap_sign(n, L):
    """Action of the input <-> output swap on the (n, n) spin-L basis: |j1 j2 L> -> (-1)^{j1+j2-L} |j2 j1 L>."""
    B = basis(n, n, L)
    idx = {b: i for i, b in enumerate(B)}
    S = np.zeros((len(B), len(B)))
    for c, (j1, j2) in enumerate(B):
        S[idx[(j2, j1)], c] = (-1) ** ((j1 + j2 - L) % 2)
    return S


# ---------------- Johnson-scheme (Hahn) data, exact ----------------
@lru_cache(maxsize=None)
def omega(m, j, s):
    """Spherical function of J(2m, m): eigenvalue of the compressed single-arrangement word with s swaps (F) and
    m - s identities on the GL(3)-isotypic component V_{(2m-j, j)} of Sym^m (x) Sym^m."""
    tot = Fraction(0)
    for h in range(0, j + 1):
        tot += Fraction((-1) ** h * math.comb(j, h) * math.comb(2 * m + 1 - j, h) * math.comb(s, h), math.comb(m, h) ** 2)
    return tot


@lru_cache(maxsize=None)
def mu(m, j):
    return math.comb(2 * m, j) - (math.comb(2 * m, j - 1) if j >= 1 else 0)
