"""Word-span decomposition SDP on Sym^n (x) Sym^n at a = 3/4.

Scaled problem (everything times 4^n):  K' = sum_{k+s+t=n} C(n;k,s,t) p^k (p-q)^s (p-r)^t   [symbol of 4^n K_n],
I = p^n.  Decomposition  K' - tau I = Q + Gamma(Q),  Gamma: q <-> r.
Parametrisation of Q in the word span (words with letters I, A = (I-F)/2, B = (I-Omega)/2):
   Q = sum_{s>t} C [ y_{st} W(s,t) + (1 - y_{st}) W(t,s) ] + 1/2 sum_{s=t>=1} C W(s,s) + (1 - tau)/2 I,
W(s,t) = symbol p^{n-s-t} (p-q)^s (p-r)^t  (sum over arrangements absorbed in C = C(n; n-s-t, s, t)).
Greedy certificate: y = 1, tau = 1.
"""
import math
from fractions import Fraction
import numpy as np
from cpl import basis, omega
from words import level_data


def multinom(n, k, s, t):
    return math.factorial(n) // (math.factorial(k) * math.factorial(s) * math.factorial(t))


def poly_W(n, s, t):
    """Monomial coefficients {(a,b,c): coeff} of p^{n-s-t} (p-q)^s (p-r)^t."""
    k = n - s - t
    out = {}
    for b in range(s + 1):
        for c in range(t + 1):
            key = (k + (s - b) + (t - c), b, c)
            out[key] = out.get(key, 0) + math.comb(s, b) * math.comb(t, c) * (-1) ** (b + c)
    return out


class Blocks:
    """Per spin L: level data; turns monomial polynomials into blocks."""

    def __init__(self, n, Ls=None):
        self.n = n
        self.Ls = list(range(0, 2 * n + 1)) if Ls is None else Ls
        self.ld = {L: level_data(n, L) for L in self.Ls}
        self.om = {}
        for m in range(0, n + 1):
            self.om[m] = np.array([[float(omega(m, j, s)) for s in range(m + 1)] for j in range(m + 1)])

    def gamma_coeffs(self, poly):
        """poly {(a,b,c): coeff}  ->  {c: array over j of sum_b coeff(a,b,c) omega_j^{(n-c)}(b)}"""
        n = self.n
        g = {}
        for (a, b, c), v in poly.items():
            m = n - c
            if c not in g:
                g[c] = np.zeros(m + 1)
            g[c] += float(v) * self.om[m][:, b]
        return g

    def block(self, poly, L):
        g = self.gamma_coeffs(poly)
        d = len(basis(self.n, self.n, L))
        M = np.zeros((d, d))
        for c, gam in g.items():
            Y, lab = self.ld[L][c]
            if Y.shape[0] == 0:
                continue
            M += (Y.T * gam[lab]) @ Y
        return (M + M.T) / 2


def pairs(n):
    return [(s, t) for t in range(0, n + 1) for s in range(t + 1, n + 1 - t)]


def add(P, Q, c=1):
    out = dict(P)
    for k, v in Q.items():
        out[k] = out.get(k, 0) + c * v
    return out


def Q0_poly(n):
    """Fixed part of Q at y = 0, without the identity part:  sum_{s>t} C W(t,s) + 1/2 sum_{s=t>=1} C W(s,s)."""
    P = {}
    for (s, t) in pairs(n):
        P = add(P, poly_W(n, t, s), multinom(n, n - s - t, s, t))
    for s in range(1, n // 2 + 1):
        P = add(P, poly_W(n, s, s), Fraction(multinom(n, n - 2 * s, s, s), 2))
    return P


def F_poly(n, s, t):
    C = multinom(n, n - s - t, s, t)
    return add({k: C * v for k, v in poly_W(n, s, t).items()}, poly_W(n, t, s), -C)


def K_poly(n):
    P = {}
    for s in range(n + 1):
        for t in range(n + 1 - s):
            P = add(P, poly_W(n, s, t), multinom(n, n - s - t, s, t))
    return P


def ident_poly(n):
    return {(n, 0, 0): 1}
