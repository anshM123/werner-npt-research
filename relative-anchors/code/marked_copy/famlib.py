"""The orbit calculus of mlib.py with the size parameters symbolic: q0 = j - i, t = u - j - i (u = k - 1, so the
number of unpaired copies of U is m = q0 + t) and i, each fixed or symbolic, and sigma always symbolic.  Entries are
rational functions in flintfield.Field (python-flint fmpq_mpoly).

Formal continuation (as in orbit3.Sym3 / symf.SymF of relative-anchors/code): the Dicke weights are the rational
functions C(m - r, q0 + d) / C(m, q0 + sh) (times 2^{-pD}), the raising and lowering factors (q0 + d + 1) and
(t - r - d + 1) are polynomials, and range checks are kept only for fixed parameters.  A key outside the range of a
symbolic parameter ("phantom") has weight 0, and the factors that would take it back into the range vanish, so the
part of a vector inside the range evolves exactly as under the true operators and phantoms never contribute.  E^-inverse uses the nodes i' = i + e of the window e = -2..2 that are genuine for the fixed parameters
(for symbolic q0 >= 4, t >= 4, i >= 2 every window node is a genuine eigenvalue of both parts, the 2 x 2 formulas apply
at every node, and distinct nodes have distinct values).
The numeric form of mlib.py at an integer point of the family equals N0 times the value of this one there,
N0 = 2^i C(m, q0 + sh) > 0 (norm0)."""
import itertools
from fractions import Fraction as Fr
from math import comb
import flint
from flintfield import Field


def KF(K, x):
    """a rational constant (int or Fraction) as an element of the field K"""
    x = Fr(x)
    return K(flint.fmpq(x.numerator, x.denominator))

NAMES = ['O0', 'O1', 'A00', 'A01', 'A10', 'A11',
         'B000', 'B001', 'B010', 'B011', 'B100', 'B101', 'B110', 'B111', 'G0', 'G1']
SMIN = {'q0': 4, 't': 4, 'i': 2}


class FamCalc:
    def __init__(self, q0=None, t=None, i=None, bscale=1):
        names = [nm for nm, v in (('q0', q0), ('t', t), ('i', i)) if v is None] + ['s']
        self.syms = names
        self.K = Field(names)
        g = dict(zip(names, self.K.gens))
        self.fq0, self.ft, self.fi = q0, t, i
        self.q0 = g['q0'] if q0 is None else self.K(q0)
        self.t = g['t'] if t is None else self.K(t)
        self.i = g['i'] if i is None else self.K(i)
        self.sig = g['s']
        self.one = self.K(1)
        self.zero = self.K(0)
        self.m = self.q0 + self.t
        self.u = self.m + 2 * self.i
        self.bscale = KF(self.K, bscale)
        self.sym = True
        # normalisation shift: C(m, q0 + sh) must be nonzero on the family
        sh = 0
        if q0 is not None and q0 < 0:
            sh = -q0
        if t is not None and t < 0:
            sh = t
        self.sh = sh
        self._einv = {}

    # ---- conditions that are decided by the fixed parameters (symbolic ones are at least SMIN)
    def minval(self, nm):
        v = {'q0': self.fq0, 't': self.ft, 'i': self.fi}[nm]
        return SMIN[nm] if v is None else v

    def ige(self, n):
        return self.fi is None or self.fi >= n

    def mge(self, n):
        if self.fq0 is None or self.ft is None:
            return True
        return self.fq0 + self.ft >= n

    def inrange(self, d, r):
        if self.fq0 is not None and self.fq0 + d < 0:
            return False
        if self.ft is not None and d > self.ft - r:
            return False
        return True

    def upf(self, d):
        return self.q0 + (d + 1)

    def dnf(self, d, r):
        return self.t - (r + d - 1)

    def weight(self, pD, r, d):
        K = self.K
        m = self.m
        q0 = self.q0 + self.sh
        d = d - self.sh
        w = KF(K, Fr(1, 2 ** pD))
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

    # ---- E^-inverse data at level j + delta
    def einv_data(self, delta):
        if delta in self._einv:
            return self._einv[delta]
        K, s = self.K, self.sig
        q0, t, i = self.q0, self.t, self.i
        half = KF(K, Fr(1, 2))

        def genuine(e, part):
            if self.fi is not None and self.fi + e < 0:
                return False
            off = 0 if part == 'A' else 1
            if self.fq0 is not None and e > self.fq0 + delta - off:
                return False
            if self.ft is not None and e > self.ft - delta + off:
                return False
            return True

        nA = [e for e in range(-2, 3) if genuine(e, 'A')]
        nB = [e for e in range(-2, 3) if genuine(e, 'B')]
        vals = {}
        for e in sorted(set(nA) | set(nB)):
            mu = (t - (delta - 1 + e)) * (q0 + (delta - e))
            nu = (t - (delta - 2 + e)) * (q0 + (delta - 1 - e))
            umc = t + i - delta                      # u - c
            a = 2 * s + 2 * umc + half * mu
            ee = 2 * (umc + 1) + half * s + half * nu
            D = a * ee - s * mu * KF(K, Fr(1, 4))
            kind = 'both' if (e in nA and e in nB) else ('A' if e in nA else 'B')
            vals[e] = (kind, a, ee, D)

        def fval(kind_f, e):
            kind, a, ee, D = vals[e]
            if kind_f == 1:
                return ee / D if kind == 'both' else self.one / a
            if kind_f == 2:
                return (-half * s) / D if kind == 'both' else self.zero
            return (s * a) / D if kind == 'both' else s / ee

        def lam(e, part):
            # DU on U-level ell = j + delta (A) or j + delta - 1 (B): (u - ell - i')(ell + 1 - i')
            off = 0 if part == 'A' else 1
            return (t - (delta - off + e)) * (q0 + (delta - off + 1 - e))

        def newton(kind_f, nodes, part):
            lamv = [lam(e, part) for e in nodes]
            coefs, diag = [], []
            for n, e in enumerate(nodes):
                new = [fval(kind_f, e)]
                for l in range(1, n + 1):
                    new.append((new[l - 1] - diag[l - 1]) / (lamv[n] - lamv[n - l]))
                diag = new
                coefs.append(new[n])
            return lamv, coefs

        data = {'A1': newton(1, nA, 'A'), 'B2': newton(2, nB, 'B'), 'B3': newton(3, nB, 'B')}
        self._einv[delta] = data
        return data


class GFrame:
    def __init__(self, calc, labels):
        self.C = calc
        self.labels = labels
        self.pairs = sorted({l[1] for l in labels if l[0] == 'P'})
        self.rs = [l for l in labels if l[0] == 'R']
        self.r = len(self.rs)
        self.pD = len(self.pairs)
        self.pos = {l: n for n, l in enumerate(labels)}
        self.upos = [n for n, l in enumerate(labels) if l != 'M']

    def inrange(self, d):
        return self.C.inrange(d, self.r)

    def delta(self, key):
        L, d = key
        return sum(L) + d - self.pD

    @staticmethod
    def addto(acc, key, val):
        nv = acc.get(key)
        nv = val if nv is None else nv + val
        if nv == 0:
            acc.pop(key, None)
        else:
            acc[key] = nv

    def add(self, acc, v, c=None):
        for key, val in v.items():
            self.addto(acc, key, val if c is None else c * val)
        return acc

    def scale(self, v, c):
        out = {}
        for key, val in v.items():
            x = c * val
            if x != 0:
                out[key] = x
        return out

    def jUp(self, v):
        out = {}
        for (L, d), c in v.items():
            for s in self.upos:
                if L[s] == 0:
                    self.addto(out, (L[:s] + (1,) + L[s + 1:], d), c)
            if self.inrange(d + 1):
                self.addto(out, (L, d + 1), c * self.C.upf(d))
        return out

    def jUm(self, v):
        out = {}
        for (L, d), c in v.items():
            for s in self.upos:
                if L[s] == 1:
                    self.addto(out, (L[:s] + (0,) + L[s + 1:], d), c)
            if self.inrange(d - 1):
                self.addto(out, (L, d - 1), c * self.C.dnf(d, self.r))
        return out

    def flipM(self, v, frm, to):
        return {(((to,) + L[1:]), d): c for (L, d), c in v.items() if L[0] == frm}

    def jp(self, v):
        return self.add(self.jUp(v), self.flipM(v, 0, 1))

    def jm(self, v):
        return self.add(self.jUm(v), self.flipM(v, 1, 0))

    def sm(self, lab, v):
        s = self.pos[lab]
        return {(L[:s] + (0,) + L[s + 1:], d): c for (L, d), c in v.items() if L[s] == 1}

    def nop(self, lab, v):
        s = self.pos[lab]
        return {key: c for key, c in v.items() if key[0][s] == 1}

    def Nop(self, lab, v):
        s = self.pos[lab]
        h, th = KF(self.C.K, Fr(1, 2)), KF(self.C.K, Fr(3, 2))
        return {(L, d): c * (h if L[s] == 1 else th) for (L, d), c in v.items()}

    def t(self, lab, v):
        out = self.scale(self.sm(lab, self.jp(v)), KF(self.C.K, Fr(1, 2)))
        return self.add(out, self.Nop(lab, v))

    def Pi(self, v):
        s = self.C.sig
        return {(L, d): (c * s if L[0] == 1 else c) for (L, d), c in v.items()}

    def DU(self, v):
        return self.jUm(self.jUp(v))

    def newton(self, v, lamv, coefs):
        out = {}
        w = dict(v)
        for n in range(len(lamv)):
            if not w:
                break
            self.add(out, w, coefs[n])
            nw = self.DU(w)
            self.add(nw, w, -lamv[n])
            w = nw
        return out

    def einv(self, v):
        if not v:
            return {}
        lv = {self.delta(key) for key in v}
        assert len(lv) == 1, lv
        data = self.C.einv_data(lv.pop())
        vA = {key: x for key, x in v.items() if key[0][0] == 0}
        vB = {key: x for key, x in v.items() if key[0][0] == 1}
        out = {}
        if vA:
            self.add(out, self.newton(vA, *data['A1']))
            jv = self.jUm(vA)
            if jv:
                self.add(out, self.flipM(self.newton(jv, *data['B2']), 0, 1))
        if vB:
            self.add(out, self.newton(vB, *data['B3']))
            y = self.newton(vB, *data['B2'])
            if y:
                self.add(out, self.jUp(self.flipM(y, 1, 0)))
        return out

    def dot(self, a, b):
        if len(a) > len(b):
            a, b = b, a
        s = self.C.zero
        for key, x in a.items():
            y = b.get(key)
            if y is not None:
                s = s + x * y * self.C.weight(self.pD, self.r, key[1])
        return s

    def _expand(self, fixed, qtot):
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
                    L[self.pos['P%sa' % p]] = 1 - sbit
                    L[self.pos['P%sb' % p]] = sbit
                    coef *= (1 if sbit == 0 else -1)
                for r, x in zip(free_rs, rl):
                    L[self.pos[r]] = x
                d = qtot - sum(rl)
                if not self.inrange(d):
                    continue
                self.addto(out, (tuple(L), d), self.C.one * coef)
        return out

    def psi(self, name, s):
        if name[0] == 'O':
            if s != 'M':
                return {}
            y0 = int(name[1])
            return self._expand({'M': y0}, -y0)
        if name[0] == 'A':
            if s[0] != 'R':
                return {}
            x, y0 = int(name[1]), int(name[2])
            return self._expand({s: x, 'M': y0}, -x - y0)
        if name[0] == 'B':
            if s[0] != 'P':
                return {}
            x, y, y0 = int(name[1]), int(name[2]), int(name[3])
            p = s[1]
            first = s[2] == 'a'
            partner = 'P%s%s' % (p, 'b' if first else 'a')
            v = self._expand({s: x, partner: y, 'M': y0}, 1 - x - y - y0)
            return v if first else self.scale(v, -1)
        return {}

    def phi(self, name):
        if name[0] != 'G':
            return {}
        z0 = int(name[1])
        return self._expand({'M': z0}, 1 - z0)


def frame_labels(*sites):
    labs = ['M']
    for s in sites:
        if s == 'M':
            continue
        if s[0] == 'P':
            for t in ('P%sa' % s[1], 'P%sb' % s[1]):
                if t not in labs:
                    labs.append(t)
        elif s not in labs:
            labs.append(s)
    return labs


def orbits(C):
    m, i = C.m, C.i
    sr = [('M', 1)]
    if C.mge(1):
        sr.append(('R0', m))
    if C.ige(1):
        sr.append(('P0a', 2 * i))
    pr = [('M', 'M', 1)]
    if C.mge(1):
        pr += [('M', 'R0', m), ('R0', 'M', m), ('R0', 'R0', m)]
    if C.mge(2):
        pr.append(('R0', 'R1', m * (m - 1)))
    if C.ige(1):
        pr += [('M', 'P0a', 2 * i), ('P0a', 'M', 2 * i), ('P0a', 'P0a', 2 * i), ('P0a', 'P0b', 2 * i)]
        if C.mge(1):
            pr += [('R0', 'P0a', 2 * i * m), ('P0a', 'R0', 2 * i * m)]
    if C.ige(2):
        pr.append(('P0a', 'P1a', 2 * i * (2 * i - 2)))
    return sr, pr


def present_names(C):
    """the twisted vectors that are nonzero on the family (decided at its minimal point; presence is monotone)"""
    q0, t, i = C.minval('q0'), C.minval('t'), C.minval('i')
    m = q0 + t
    u = m + 2 * i
    k, j = u + 1, q0 + i
    out = []
    for nm in NAMES:
        if nm[0] == 'O':
            y0 = int(nm[1]); ok = 0 <= q0 - y0 <= m
        elif nm[0] == 'A':
            x, y0 = int(nm[1]), int(nm[2]); ok = m >= 1 and 0 <= q0 - y0 - x <= m - 1
        elif nm[0] == 'B':
            x, y, y0 = int(nm[1]), int(nm[2]), int(nm[3]); ok = i >= 1 and 0 <= q0 + 1 - x - y - y0 <= m
        else:
            z0 = int(nm[1]); ok = j + 1 <= k and 0 <= q0 + 1 - z0 <= m
        if ok:
            out.append(nm)
    return out


def reduced(C):
    names = present_names(C)
    n = len(NAMES)
    sr, pr = orbits(C)
    sig = C.sig
    Q = [[C.zero] * n for _ in range(n)]
    b = C.bscale
    s_of = lambda lab: sig if lab == 'M' else C.one
    half = KF(C.K, Fr(1, 2))
    for (a, sz) in sr:
        fr = GFrame(C, frame_labels(a))
        comps = {nm: fr.psi(nm, a) for nm in names}
        phis = {nm: fr.phi(nm) for nm in names}
        sa = s_of(a)
        Cv = {}
        for nm in names:
            ph = phis[nm]
            if ph:
                v = fr.scale(fr.sm(a, fr.Pi(ph)), -1)
                fr.add(v, fr.scale(fr.jm(fr.einv(fr.t(a, ph))), sa * b))
                Cv[nm] = v
        for na in names:
            pa = comps[na]
            if not pa:
                continue
            x = NAMES.index(na)
            Pipa = fr.Pi(pa)
            for nb in names:
                y = NAMES.index(nb)
                pb = comps[nb]
                if pb:
                    Q[x][y] += sz * sa * (2 * fr.dot(Pipa, pb) - fr.dot(Pipa, fr.nop(a, pb)))
                if nb in Cv:
                    val = sz * fr.dot(pa, Cv[nb])
                    Q[x][y] += val
                    Q[y][x] += val
    for (a, bb, sz) in pr:
        fr = GFrame(C, frame_labels(a, bb))
        compa = {nm: fr.psi(nm, a) for nm in names}
        compb = {nm: fr.psi(nm, bb) for nm in names}
        sab = s_of(a) * s_of(bb)
        coef = half - (1 if a == bb else 0)
        Ec, smA = {}, {}
        for nb in names:
            pb = compb[nb]
            if pb:
                Ec[nb] = fr.einv(fr.t(a, pb))
                smA[nb] = fr.Pi(fr.sm(a, pb))
        for na in names:
            pa = compa[na]
            if not pa:
                continue
            x = NAMES.index(na)
            Tb = fr.t(bb, pa)
            smb = fr.sm(bb, pa)
            for nb in Ec:
                y = NAMES.index(nb)
                Q[x][y] += sz * (coef * sab * fr.dot(smb, smA[nb]) - b * sab * fr.dot(Tb, Ec[nb]))
    fr = GFrame(C, ['M'])
    ph = {nm: fr.phi(nm) for nm in names}
    for na in names:
        for nb in names:
            if ph[na] and ph[nb]:
                x, y = NAMES.index(na), NAMES.index(nb)
                Q[x][y] += 2 * fr.dot(fr.Pi(ph[na]), ph[nb])
                Ua, Ub = fr.jp(ph[na]), fr.jp(ph[nb])
                if Ua and Ub:
                    Q[x][y] += -b * fr.dot(Ua, fr.einv(Ub))
    idx = [NAMES.index(nm) for nm in names]
    return names, [[Q[x][y] for y in idx] for x in idx]


def norm0(C, q0, t, i):
    """N0 = 2^i C(m, q0 + sh) at an integer point: the numeric form of mlib equals N0 times this one"""
    m = q0 + t
    return (2 ** i) * comb(m, q0 + C.sh)
