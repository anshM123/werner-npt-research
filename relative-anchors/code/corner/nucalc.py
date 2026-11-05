"""Backends and helpers for the corner piece (the (u - i, i) part of the sym case I letter-2 forms).

Variables.  The generic-family calculus orbit2.SymF2 works over QQ(q0, t, i, tau).  Every frame quantity depends on
(i, tau) only through nu = 2i + tau (the constant of calN, and hence every Newton denominator), while i and tau also
enter the H-orbit sizes, which are polynomials of degree <= 2 in (i, tau).  SymNu is the same calculus over
QQ(q0, t, nu) with tau a FIXED (formal) rational number and i = (nu - tau)/2.  All orbit classes are kept for every
value of tau (their sizes vanish where the class is empty), so the result is the specialisation at that tau of one
rational function of (q0, t, nu, tau), polynomial of degree <= 2 in tau with tau-free denominators.

Q7 = the (u - i, i) form Q8 in the basis X1, X2, X3, X4, G0, D0, D1 (C1 dropped; kappa = X1 + X2 + C1 spans the
kernel of Q8, so Q8 is congruent to Q7 (+) 0)."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
from orbit2 import SymF2, reduced2

V8c = {'X1': {'A1': 1, 'B10': 1}, 'X2': {'A0': 1, 'B01': -1}, 'X3': {'A1': 1, 'B11': -1}, 'X4': {'A0': 1, 'B00': 1},
       'G0': {'G0': 1}, 'C1': {'C1': 1}, 'D0': {'D0': 1}, 'D1': {'D1': 1}}
Q7NAMES = ['X1', 'X2', 'X3', 'X4', 'G0', 'D0', 'D1']


class SymNu(SymF2):
    """QQ(q0, t, nu) backend, tau fixed (any rational, formal), i = (nu - tau)/2; generic family (q0 >= 2, t >= 4,
    i >= 2: node offsets -2..2)."""

    def __init__(self, tau):
        from flintfield import Field
        self.syms = ['q0', 't', 'nu']
        self.Kf = Field(self.syms)
        q0, t, nu = self.Kf.gens
        self.fq0 = self.ft = self.fi = self.ftau = None
        self.q0, self.t = q0, t
        self.tau = self.Kf(tau)
        self.i = (nu - tau) / 2
        self.nu = nu
        self.m = self.q0 + self.t
        self.u = self.m + 2 * self.i
        self.sh = 0
        self.seen = set()

    def tge(self, n):
        return True        # keep the T-T pair class for every formal tau (size tau(tau-1))


def form_in_basis(be, names_out):
    names, G, Q = reduced2(be, 'sym')
    K = be.K
    idx = {nm: names.index(nm) for nm in names}

    def vec(nm):
        v = [K(0)] * len(names)
        for x, c in V8c[nm].items():
            v[idx[x]] = K(c)
        return v
    n = len(names)

    def bil(M, u, v):
        return sum((u[x] * sum((M[x][y] * v[y] for y in range(n) if v[y] != 0), K(0)) for x in range(n) if u[x] != 0), K(0))
    vs = [vec(nm) for nm in names_out]
    Qo = [[bil(Q, a, b) for b in vs] for a in vs]
    Go = [[bil(G, a, b) for b in vs] for a in vs]
    return Qo, Go


def det_field(M, K, order=None):
    """determinant by symmetric Gaussian elimination over the field (pivots must be nonzero)"""
    n = len(M)
    M = [row[:] for row in M]
    order = list(range(n)) if order is None else list(order)
    act = list(order)
    det = K(1)
    for x in order:
        piv = M[x][x]
        if piv == 0:
            raise ZeroDivisionError('zero pivot')
        det = det * piv
        act.remove(x)
        for r in act:
            f = M[r][x] / piv
            if f != 0:
                for s in act:
                    M[r][s] = M[r][s] - f * M[x][s]
    return det
