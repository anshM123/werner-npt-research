"""Orbit calculus for the equal-weight letter-1 block of (R_inf)_j: no enumeration of strings.

A vector in a frame Delta (at most 4 distinguished sites: up to two pair-sites with partners, up to two R-sites)
is a dict (letters on Delta, q) -> coefficient, meaning
    sum coef * |letters>_Delta (x) [singlets on the i - p_Delta undistinguished pairs] (x) Dicke_{R - Delta}(q).
U and D kill singlets, so  U: flip a 0 in Delta, or Dicke q -> q+1 with factor (q+1);
D: flip a 1 in Delta, or Dicke q -> q-1 with factor (n' - q + 1), n' = |R - Delta| = m - |R cap Delta|.
<(L,q),(L,q)> = 2^(i - p_Delta) C(n', q).  E_c^{-1} = sum_n f[lam_0..lam_n] prod_{l<n} (DU - lam_l),
f(lam) = 1/(alpha_c + lam/2), f[lam_0..lam_n] = (-1/2)^n / prod_l (alpha_c + lam_l/2), nodes lam_{i'} =
(k-c-i')(c+1-i') over the spin components i' that can occur.
Numeric mode: exact Fractions at integer k (any k, cost independent of k).  Symbolic mode: sympy, k symbolic.
"""
from fractions import Fraction
from math import comb
import itertools


class Num:
    """exact numeric backend at integer k"""

    def __init__(self, k, j, i):
        self.k, self.j, self.i = k, j, i
        self.m = k - 2 * i

    def C(self, n, q):
        if q < 0 or n < 0 or q > n:
            return 0
        return comb(n, q)

    def nodes(self, c):
        k, i = self.k, self.i
        return [ip for ip in range(0, min(c, k - c) + 1)]

    def lam(self, c, ip):
        return Fraction((self.k - c - ip) * (c + 1 - ip))

    def alpha(self, c):
        return Fraction(3 * self.k - 2 * c, 2)

    def F(self, x):
        return Fraction(x)

    def R(self, a, b):
        return Fraction(a, b)

    def le(self, a, b):
        return a <= b


class Sym:
    """symbolic backend: k a symbol (rational function field QQ(k)), j and i fixed integers; valid for k large
    enough that every orbit type exists and min(c, k - c) = c for the levels c <= j + 2 used (k >= 2j + 4)."""

    def __init__(self, j, i):
        from sympy import QQ
        from sympy.polys.fields import field
        self.K, self.kk = field('k', QQ)
        self.k, self.j, self.i = self.kk, j, i
        self.m = self.kk - 2 * i

    def C(self, n, q):
        if q < 0:
            return 0
        v = self.K(1)
        for t in range(q):
            v = v * (n - t) / (t + 1)
        return v

    def nodes(self, c):
        return list(range(0, c + 1))

    def lam(self, c, ip):
        return (self.k - c - ip) * (c + 1 - ip)

    def alpha(self, c):
        return (3 * self.k - 2 * c) / self.K(2)

    def R(self, a, b):
        return self.K(a) / b

    def le(self, a, b):
        if isinstance(b, int):
            return a <= b
        return True


class Frame:
    def __init__(self, be, labels):
        """labels: list like ['P0a','P0b','R0']; pair labels must come with partners"""
        self.be = be
        self.labels = labels
        self.pairs = sorted({l[1] for l in labels if l[0] == 'P'})
        self.rs = [l for l in labels if l[0] == 'R']
        self.pD = len(self.pairs)
        self.nprime = be.m - len(self.rs)
        self.undist = be.i - self.pD
        self.pos = {l: n for n, l in enumerate(labels)}

    def level(self, L, q):
        return sum(L) + q + self.undist

    # ---- vector ops
    def add(self, acc, v, c=1):
        for key, val in v.items():
            nv = acc.get(key, 0) + c * val
            acc[key] = nv
        return acc

    def clean(self, v):
        return {k_: x for k_, x in v.items() if x != 0}

    def U(self, v):
        out = {}
        for (L, q), c in v.items():
            for s in range(len(L)):
                if L[s] == 0:
                    L2 = L[:s] + (1,) + L[s + 1:]
                    out[(L2, q)] = out.get((L2, q), 0) + c
            if self.be.le(q + 1, self.nprime):
                out[(L, q + 1)] = out.get((L, q + 1), 0) + c * (q + 1)
        return self.clean(out)

    def D(self, v):
        out = {}
        for (L, q), c in v.items():
            for s in range(len(L)):
                if L[s] == 1:
                    L2 = L[:s] + (0,) + L[s + 1:]
                    out[(L2, q)] = out.get((L2, q), 0) + c
            if q >= 1:
                out[(L, q - 1)] = out.get((L, q - 1), 0) + c * (self.nprime - q + 1)
        return self.clean(out)

    def sm(self, lab, v):
        s = self.pos[lab]
        return {(L[:s] + (0,) + L[s + 1:], q): c for (L, q), c in v.items() if L[s] == 1}

    def nop(self, lab, v):
        s = self.pos[lab]
        return {(L, q): c for (L, q), c in v.items() if L[s] == 1}

    def Nop(self, lab, v):
        s = self.pos[lab]
        return {(L, q): c * (self.be.R(1, 2) if L[s] == 1 else self.be.R(3, 2)) for (L, q), c in v.items()}

    def T(self, lab, v):
        out = {key: self.be.R(1, 2) * c for key, c in self.sm(lab, self.U(v)).items()}
        return self.clean(self.add(out, self.Nop(lab, v)))

    def levelof(self, v):
        lv = {self.level(L, q) for (L, q) in v}
        assert len(lv) <= 1, lv
        return lv.pop() if lv else None

    def Einv(self, v):
        if not v:
            return {}
        c = self.levelof(v)
        be = self.be
        nodes = be.nodes(c)
        lam = [be.lam(c, ip) for ip in nodes]
        al = be.alpha(c)
        out = {}
        w = dict(v)
        den = 1
        for n in range(len(lam)):
            den = den * (al + lam[n] / 2)
            coef = self.be.R(-1, 2) ** n / den
            self.add(out, w, coef)
            if n < len(lam) - 1:
                DUw = self.D(self.U(w))
                w = self.clean(self.add(DUw, w, -lam[n]))
        return self.clean(out)

    def dot(self, u, v):
        s = 0
        for key, c in u.items():
            d = v.get(key)
            if d is not None:
                L, q = key
                s += c * d * (2 ** self.undist) * self.be.C(self.nprime, q)
        return s

    # ---- the basis vectors in this frame
    def _expand(self, fixed, qtot):
        """fixed: dict label->letter for some labels; other pair labels -> singlets; other R labels -> Dicke.
        Returns vector with total Dicke count qtot over (R cap Delta minus fixed) + R'."""
        free_pairs = [p for p in self.pairs if ('P%sa' % p) not in fixed]
        free_rs = [r for r in self.rs if r not in fixed]
        out = {}
        for sing in itertools.product((0, 1), repeat=len(free_pairs)):
            for rl in itertools.product((0, 1), repeat=len(free_rs)):
                L = [None] * len(self.labels)
                coef = 1
                for lab, x in fixed.items():
                    L[self.pos[lab]] = x
                for p, sbit in zip(free_pairs, sing):
                    # singlet |1 at a, 0 at b> - |0 at a, 1 at b>
                    L[self.pos['P%sa' % p]] = 1 - sbit
                    L[self.pos['P%sb' % p]] = sbit
                    coef *= (1 if sbit == 0 else -1)
                for r, x in zip(free_rs, rl):
                    L[self.pos[r]] = x
                q = qtot - sum(rl)
                if q < 0 or not self.be.le(q, self.nprime):
                    continue
                key = (tuple(L), q)
                out[key] = out.get(key, 0) + coef
        return self.clean(out)

    def psi(self, name, s):
        """component at site label s of basis vector `name`"""
        be = self.be
        j, i = be.j, be.i
        if name in ('A0', 'A1'):
            x = int(name[1])
            if s[0] != 'R':
                return {}
            return self._expand({s: x}, j - i - x)
        if name[0] == 'B':
            x, y = int(name[1]), int(name[2])
            if s[0] != 'P':
                return {}
            p = s[1]
            first = s[2] == 'a'
            partner = 'P%s%s' % (p, 'b' if first else 'a')
            v = self._expand({s: x, partner: y}, j - i + 1 - x - y)
            return {key: (c if first else -c) for key, c in v.items()}
        return {}

    def phi(self, name):
        if name != 'G0':
            return {}
        be = self.be
        return self._expand({}, be.j + 1 - be.i)


NAMES = ['A0', 'A1', 'B00', 'B01', 'B10', 'B11', 'G0']


def frame_labels(*sites):
    labs = []
    for s in sites:
        if s[0] == 'P':
            for t in ('P%sa' % s[1], 'P%sb' % s[1]):
                if t not in labs:
                    labs.append(t)
        elif s not in labs:
            labs.append(s)
    return labs


def reduced_orbit(be):
    k, j, i, m = be.k, be.j, be.i, be.m
    pr = []
    big = not isinstance(m, int)
    if big or m >= 1:
        pr.append(('R0', 'R0', m))
    if big or m >= 2:
        pr.append(('R0', 'R1', m * (m - 1)))
    if (big or m >= 1) and i >= 1:
        pr.append(('R0', 'P0a', 2 * i * m))
        pr.append(('P0a', 'R0', 2 * i * m))
    if i >= 1:
        pr.append(('P0a', 'P0a', 2 * i))
        pr.append(('P0a', 'P0b', 2 * i))
    if i >= 2:
        pr.append(('P0a', 'P1a', 2 * i * (2 * i - 2)))
    sr = []
    if big or m >= 1:
        sr.append(('R0', m))
    if i >= 1:
        sr.append(('P0a', 2 * i))
    n = len(NAMES)
    G = [[0] * n for _ in range(n)]
    Q = [[0] * n for _ in range(n)]
    # site terms
    for (a, sz) in sr:
        fr = Frame(be, frame_labels(a))
        comps = {nm: fr.psi(nm, a) for nm in NAMES}
        phis = {nm: fr.phi(nm) for nm in NAMES}
        Cv = {}
        for nm in NAMES:
            ph = phis[nm]
            if ph:
                v = {key: -c for key, c in fr.sm(a, ph).items()}
                fr.add(v, fr.D(fr.Einv(fr.T(a, ph))))
                Cv[nm] = fr.clean(v)
        for x, na in enumerate(NAMES):
            for y, nb in enumerate(NAMES):
                pa, pb = comps[na], comps[nb]
                if pa and pb:
                    G[x][y] += sz * fr.dot(pa, pb)
                    Q[x][y] += sz * (2 * fr.dot(pa, pb) - fr.dot(pa, fr.nop(a, pb)))
                if pa and nb in Cv:
                    Q[x][y] += sz * fr.dot(pa, Cv[nb])
                if pb and na in Cv:
                    Q[x][y] += sz * fr.dot(pb, Cv[na])
    # pair terms
    for (a, b, sz) in pr:
        fr = Frame(be, frame_labels(a, b))
        compa = {nm: fr.psi(nm, a) for nm in NAMES}
        compb = {nm: fr.psi(nm, b) for nm in NAMES}
        Ecache = {}
        for nb in NAMES:
            pb = compb[nb]
            if pb:
                Ecache[nb] = fr.Einv(fr.T(a, pb))
        for x, na in enumerate(NAMES):
            pa = compa[na]
            if not pa:
                continue
            Tb_pa = fr.T(b, pa)
            smb_pa = fr.sm(b, pa)
            for y, nb in enumerate(NAMES):
                pb = compb[nb]
                if not pb:
                    continue
                coef = be.R(1, 2) - (1 if a == b else 0)
                Q[x][y] += sz * coef * fr.dot(smb_pa, fr.sm(a, pb))
                Q[x][y] -= sz * fr.dot(Tb_pa, Ecache[nb])
    # G0 terms
    fr = Frame(be, [])
    ph = {nm: fr.phi(nm) for nm in NAMES}
    for x, na in enumerate(NAMES):
        for y, nb in enumerate(NAMES):
            if ph[na] and ph[nb]:
                G[x][y] += fr.dot(ph[na], ph[nb])
                Q[x][y] += 2 * fr.dot(ph[na], ph[nb])
                Ua, Ub = fr.U(ph[na]), fr.U(ph[nb])
                if Ua and Ub:
                    Q[x][y] -= fr.dot(Ua, fr.Einv(Ub))
    keep = [x for x in range(n) if G[x][x] != 0]
    return [NAMES[x] for x in keep], [[G[x][y] for y in keep] for x in keep], [[Q[x][y] for y in keep] for x in keep]
