"""flint-based exact-node symbolic backend (see symt.py for the validity argument).  Parameters q0, t, i:
fixed int or None (symbolic; then q0 >= 2, t >= 4, i >= 2 are assumed, so that every Newton node is a genuine
eigenvalue).  sh: normalisation shift (base binomial C(m, q0 + sh)); sh = 1 for q0 = -1, sh = -1 for t = -1."""
import sys
from flintfield import Field


class SymF:
    def __init__(self, q0=None, t=None, i=None):
        names = [nm for nm, v in (('q0', q0), ('t', t), ('i', i)) if v is None]
        self.syms = names
        self.Kf = Field(names if names else ['z'])
        g = dict(zip(names, self.Kf.gens))
        self.fq0, self.ft, self.fi = q0, t, i
        self.q0 = g.get('q0', q0)
        self.t = g.get('t', t)
        self.i = g.get('i', i)
        self.m = self.q0 + self.t
        self.k = self.m + 2 * self.i
        self.j = self.q0 + self.i
        self.sh = 1 if q0 == -1 else (-1 if t == -1 else 0)

    def K(self, x):
        return self.Kf(x)

    def R(self, a, b):
        return self.Kf(a) / b

    def level_pD(self, sumL, d, pD):
        return ('off', sumL + d - pD)

    def nodes(self, c):
        delta = c[1]
        lo = -2 if self.fi is None else max(-2, -self.fi)
        hi = 2
        if self.fq0 is not None:
            hi = min(hi, self.fq0 + delta)
        if self.ft is not None:
            hi = min(hi, self.ft - delta)
        return list(range(lo, hi + 1))

    def lam(self, c, e):
        cc = self.j + c[1]
        ip = self.i + e
        return (self.k - cc - ip) * (cc + 1 - ip)

    def alpha(self, c):
        cc = self.j + c[1]
        return (3 * self.k - 2 * cc) / self.Kf(2) if not isinstance(cc, int) or not isinstance(self.k, int) else self.Kf(3 * self.k - 2 * cc) / 2

    def inrange(self, d, r):
        return True

    def weight(self, pD, r, d):
        m = self.m
        q0 = self.q0 + self.sh
        d = d - self.sh
        w = self.Kf(1) / 2 ** pD
        if d >= 0:
            for s in range(1, d + 1):
                w = w / (q0 + s)
        else:
            for s in range(0, -d):
                w = w * (q0 - s)
        s_ = r + d
        if s_ >= 0:
            for s in range(0, s_):
                w = w * (m - q0 - s)
        else:
            for s in range(1, -s_ + 1):
                w = w / (m - q0 + s)
        for s in range(0, r):
            w = w / (m - s)
        return w

    def ige(self, n):
        return True if self.fi is None else self.fi >= n
