"""Orbit calculus for the equal-weight CASE I letter-2 torus blocks of (R_inf)_j (no enumeration of strings).

Block: T = copies with exactly one letter 2 (|T| = tau), U = the other u = k - tau copies = i pairs + rest R (m copies).
Base Dicke count q0 = j - tau - i (the count of C1), t = m - q0.  k = q0 + t + 2i + tau, j = q0 + i + tau.
Frame: distinguished labels Delta (U-pair sites with partners 'P0a','P0b','P1a','P1b', rest sites 'R0','R1',
T-sites 'Ta','Tb' (anti pair / two distinct T-sites) and 'Tc' (a rest T-site in anti mode)).  A frame vector is
dict (letters on Delta, d) -> coef, meaning |letters>_Delta (x) singlets on the i - p_Delta undistinguished U-pairs
(x) Dicke_{R - Delta}(q0 + d) (x) letter 2 on every undistinguished T-site.
Letters 0/1 are Boolean; letter 2 is inert for U = J+ and D = J-.  <.,.> = 2^(i - p_Delta) C(|R - Delta|, q0 + d).
E = calN + 1/2 DU on a fixed 2-pattern; calN = 3/2 #0 + 1/2 #1 + #2; DU = Johnson operator of the Boolean copies
(U and the Boolean T-sites of Delta, k_B = u + b_T copies, c' ones): lambda = (k_B - c' - i')(c' + 1 - i').
Spin window: i' in [i - 2, i + 2] (at most two row copies per frame), intersected with the valid range.
"""
from fractions import Fraction
from math import comb
import itertools

HALFF = Fraction(1, 2)


# ----------------------------------------------------------------------------------------------- backends
class Num2:
    def __init__(self, k, j, tau, i):
        self.k, self.j, self.tau, self.i = k, j, tau, i
        self.u = k - tau
        self.m = self.u - 2 * i
        self.q0 = j - tau - i
        self.t = self.m - self.q0

    def R(self, a, b):
        return Fraction(a, b)

    def level(self, ones, d, pD, bT):
        cp = ones + self.q0 + d + (self.i - pD)
        kB = self.u + bT
        return (cp, kB, bT)

    def nodes(self, lev):
        cp, kB, bT = lev
        if cp < 0 or cp > kB:
            return []
        return list(range(0, min(cp, kB - cp) + 1))

    def lam(self, lev, ip):
        cp, kB, bT = lev
        return Fraction((kB - cp - ip) * (cp + 1 - ip))

    def alpha(self, lev):
        cp, kB, bT = lev
        return Fraction(3, 2) * (kB - cp) + HALFF * cp + (self.tau - bT)

    def inrange(self, d, r):
        q = self.q0 + d
        return 0 <= q <= self.m - r

    def weight(self, pD, r, d):
        q = self.q0 + d
        n = self.m - r
        if q < 0 or q > n or n < 0:
            return 0
        return 2 ** (self.i - pD) * comb(n, q)

    def ige(self, n):
        return self.i >= n

    def tge(self, n):
        return self.tau >= n

    def mge(self, n):
        return self.m >= n


class SymF2:
    """flint QQ(params) backend; each of q0, t, i, tau fixed int or None (symbolic).  Records the level offsets met
    so that the caller can check that the node window is exact on the family."""

    def __init__(self, q0=None, t=None, i=None, tau=None):
        from flintfield import Field
        names = [nm for nm, v in (('q0', q0), ('t', t), ('i', i), ('tau', tau)) if v is None]
        self.syms = names
        self.Kf = Field(names if names else ['z'])
        g = dict(zip(names, self.Kf.gens))
        self.fq0, self.ft, self.fi, self.ftau = q0, t, i, tau
        self.q0 = g.get('q0', q0)
        self.t = g.get('t', t)
        self.i = g.get('i', i)
        self.tau = g.get('tau', tau)
        self.m = self.q0 + self.t
        self.u = self.m + 2 * self.i
        self.sh = 1 if q0 == -1 else (-1 if t == -1 else 0)
        self.seen = set()

    def K(self, x):
        return self.Kf(x)

    def R(self, a, b):
        return self.Kf(a) / b

    def level(self, ones, d, pD, bT):
        delta = ones + d - pD            # c' - i = q0 + delta
        deltap = bT - ones - d + pD      # k_B - c' - i = t + deltap
        self.seen.add((delta, deltap))
        return (delta, deltap, bT)

    def nodes(self, lev):
        delta, deltap, bT = lev
        lo = -2 if self.fi is None else max(-2, -self.fi)
        hi = 2
        if self.fq0 is not None:
            hi = min(hi, self.fq0 + delta)
        if self.ft is not None:
            hi = min(hi, self.ft + deltap)
        return list(range(lo, hi + 1))

    def _cp_kB(self, lev):
        delta, deltap, bT = lev
        return self.i + self.q0 + delta, self.u + bT

    def lam(self, lev, e):
        cp, kB = self._cp_kB(lev)
        ip = self.i + e
        return self.Kf(0) + (kB - cp - ip) * (cp + 1 - ip)

    def alpha(self, lev):
        cp, kB = self._cp_kB(lev)
        bT = lev[2]
        return (kB - cp) * self.R(3, 2) + cp * self.R(1, 2) + (self.tau - bT)

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

    def tge(self, n):
        return True if self.ftau is None else self.ftau >= n

    def mge(self, n):
        if self.fq0 is not None and self.ft is not None:
            return self.fq0 + self.ft >= n
        return True


# ----------------------------------------------------------------------------------------------- frames
class Frame2:
    def __init__(self, be, labels):
        self.be = be
        self.labels = labels
        self.pos = {l: n for n, l in enumerate(labels)}
        self.pairs = sorted({l[1] for l in labels if l[0] == 'P'})
        self.rs = [l for l in labels if l[0] == 'R']
        self.ts = [l for l in labels if l[0] == 'T']
        self.pD = len(self.pairs)
        self.nprime = be.m - len(self.rs)

    def clean(self, v):
        return {k_: x for k_, x in v.items() if x != 0}

    def add(self, acc, v, c=1):
        for key, val in v.items():
            acc[key] = acc.get(key, 0) + c * val
        return acc

    def U(self, v):
        out = {}
        for (L, d), c in v.items():
            for s in range(len(L)):
                if L[s] == 0:
                    key = (L[:s] + (1,) + L[s + 1:], d)
                    out[key] = out.get(key, 0) + c
            if self.be.inrange(d + 1, len(self.rs)):
                out[(L, d + 1)] = out.get((L, d + 1), 0) + c * (self.be.q0 + d + 1)
        return self.clean(out)

    def D(self, v):
        out = {}
        for (L, d), c in v.items():
            for s in range(len(L)):
                if L[s] == 1:
                    key = (L[:s] + (0,) + L[s + 1:], d)
                    out[key] = out.get(key, 0) + c
            if self.be.inrange(d - 1, len(self.rs)):
                out[(L, d - 1)] = out.get((L, d - 1), 0) + c * (self.nprime - (self.be.q0 + d) + 1)
        return self.clean(out)

    def sm(self, lab, l, v):
        s = self.pos[lab]
        return {(L[:s] + (0,) + L[s + 1:], d): c for (L, d), c in v.items() if L[s] == l}

    def nop(self, lab, v):
        s = self.pos[lab]
        return {(L, d): c for (L, d), c in v.items() if L[s] != 0}

    def Nop(self, lab, v):
        s = self.pos[lab]
        f = {0: self.be.R(3, 2), 1: self.be.R(1, 2), 2: self.be.R(1, 1)}
        return {(L, d): c * f[L[s]] for (L, d), c in v.items()}

    def T(self, a, v):
        lab, l = a
        out = {key: self.be.R(1, 2) * c for key, c in self.sm(lab, l, self.U(v)).items()}
        if l == 1:
            self.add(out, self.Nop(lab, v))
        else:
            s = self.pos[lab]
            self.add(out, {(L[:s] + (1,) + L[s + 1:], d): -self.be.R(1, 2) * c for (L, d), c in v.items() if L[s] == 2})
        return self.clean(out)

    def lev(self, L, d):
        ones = sum(1 for x in L if x == 1)
        bT = sum(1 for lab in self.ts if L[self.pos[lab]] != 2)
        return self.be.level(ones, d, self.pD, bT)

    def Einv(self, v):
        groups = {}
        for (L, d), c in v.items():
            pat = tuple(lab for lab in self.ts if L[self.pos[lab]] == 2)
            lv = self.lev(L, d)
            groups.setdefault((pat, lv), {})[(L, d)] = c
        out = {}
        for (pat, lv), w0 in groups.items():
            nodes = self.be.nodes(lv)
            lam = [self.be.lam(lv, e) for e in nodes]
            al = self.be.alpha(lv)
            w = dict(w0)
            den = 1
            for n in range(len(lam)):
                den = den * (al + lam[n] / 2)
                self.add(out, w, self.be.R((-1) ** n, 2 ** n) / den)
                if n < len(lam) - 1:
                    w = self.clean(self.add(self.D(self.U(w)), w, -lam[n]))
        return self.clean(out)

    def dot(self, u, v):
        s = 0
        for key, c in u.items():
            e = v.get(key)
            if e is not None:
                s += c * e * self.be.weight(self.pD, len(self.rs), key[1])
        return s

    def expand(self, fixed, dtot, sign=1):
        """fixed: label -> letter (U labels fixed by the vector; T labels default to 2)"""
        free_pairs = [p for p in self.pairs if ('P%sa' % p) not in fixed]
        free_rs = [r for r in self.rs if r not in fixed]
        out = {}
        for sing in itertools.product((0, 1), repeat=len(free_pairs)):
            for rl in itertools.product((0, 1), repeat=len(free_rs)):
                L = [None] * len(self.labels)
                coef = sign
                for lab in self.ts:
                    L[self.pos[lab]] = 2
                for lab, x in fixed.items():
                    L[self.pos[lab]] = x
                for p, sbit in zip(free_pairs, sing):
                    L[self.pos['P%sa' % p]] = 1 - sbit
                    L[self.pos['P%sb' % p]] = sbit
                    coef *= (1 if sbit == 0 else -1)
                for r, x in zip(free_rs, rl):
                    L[self.pos[r]] = x
                d = dtot - sum(rl)
                if not self.be.inrange(d, len(self.rs)):
                    continue
                key = (tuple(L), d)
                out[key] = out.get(key, 0) + coef
        return self.clean(out)

    def psi(self, name, a):
        lab, l = a
        if name in ('A0', 'A1'):
            x = int(name[1])
            return self.expand({lab: x}, -x) if (lab[0] == 'R' and l == 1) else {}
        if name[0] == 'B':
            x, y = int(name[1]), int(name[2])
            if lab[0] != 'P' or l != 1:
                return {}
            first = lab[2] == 'a'
            partner = 'P%s%s' % (lab[1], 'b' if first else 'a')
            return self.expand({lab: x, partner: y}, 1 - x - y, 1 if first else -1)
        if name == 'C1':
            return self.expand({}, 0) if (lab[0] == 'T' and l == 1) else {}
        if name in ('D0', 'D1'):
            z = int(name[1])
            return self.expand({lab: z}, 1 - z) if (lab[0] == 'T' and l == 2) else {}
        if name == 'C1a':
            if l == 1 and lab in ('Ta', 'Tb'):
                return self.expand({}, 0, 1 if lab == 'Ta' else -1)
            return {}
        if name in ('D0a', 'D1a'):
            z = int(name[1])
            if l == 2 and lab in ('Ta', 'Tb'):
                return self.expand({lab: z}, 1 - z, 1 if lab == 'Ta' else -1)
            return {}
        return {}

    def phi(self, name):
        return self.expand({}, 1) if name == 'G0' else {}


SYM = ['A0', 'A1', 'B00', 'B01', 'B10', 'B11', 'G0', 'C1', 'D0', 'D1']
ANTI = ['C1a', 'D0a', 'D1a']


def frame_labels(*sites):
    labs = []
    for s in sites:
        if s[0] == 'P':
            for t_ in ('P%sa' % s[1], 'P%sb' % s[1]):
                if t_ not in labs:
                    labs.append(t_)
        elif s not in labs:
            labs.append(s)
    return labs


def site_classes(be, mode):
    """(label, size) for row sites, by H-orbit"""
    i, m, tau = be.i, be.m, be.tau
    out = []
    if mode == 'sym':
        if be.mge(1):
            out.append(('R0', m))
        if be.ige(1):
            out.append(('P0a', 2 * i))
        out.append(('Ta', tau))
    else:
        out.append(('Ta', 2))
    return out


def pair_classes(be, mode):
    """(label_a, label_b, size) for ordered pairs of row sites, by H-orbit"""
    i, m, tau = be.i, be.m, be.tau
    pr = []
    if mode == 'sym':
        if be.mge(1):
            pr.append(('R0', 'R0', m))
        if be.mge(2):
            pr.append(('R0', 'R1', m * (m - 1)))
        if be.mge(1) and be.ige(1):
            pr.append(('R0', 'P0a', 2 * i * m))
            pr.append(('P0a', 'R0', 2 * i * m))
        if be.ige(1):
            pr.append(('P0a', 'P0a', 2 * i))
            pr.append(('P0a', 'P0b', 2 * i))
        if be.ige(2):
            pr.append(('P0a', 'P1a', 2 * i * (2 * i - 2)))
        if be.mge(1):
            pr.append(('R0', 'Ta', m * tau))
            pr.append(('Ta', 'R0', m * tau))
        if be.ige(1):
            pr.append(('P0a', 'Ta', 2 * i * tau))
            pr.append(('Ta', 'P0a', 2 * i * tau))
        pr.append(('Ta', 'Ta', tau))
        if be.tge(2):
            pr.append(('Ta', 'Tb', tau * (tau - 1)))
    else:
        pr.append(('Ta', 'Ta', 2))
        pr.append(('Ta', 'Tb', 2))
    return pr


def reduced2(be, mode):
    names = SYM if mode == 'sym' else ANTI
    n = len(names)
    G = [[0] * n for _ in range(n)]
    Q = [[0] * n for _ in range(n)]
    letters = (1, 2)
    # site terms: (2 - n) diagonal term, Gram, coupling to G0
    for (s, sz) in site_classes(be, mode):
        fr = Frame2(be, frame_labels(s))
        ph = {nm: fr.phi(nm) for nm in names}
        for l in letters:
            a = (s, l)
            comp = {nm: fr.psi(nm, a) for nm in names}
            Cv = {}
            for nm in names:
                if ph[nm]:
                    v = {key: -c for key, c in fr.sm(s, l, ph[nm]).items()}
                    fr.add(v, fr.D(fr.Einv(fr.T(a, ph[nm]))))
                    Cv[nm] = fr.clean(v)
            for x, na in enumerate(names):
                for y, nb in enumerate(names):
                    pa, pb = comp[na], comp[nb]
                    if pa and pb:
                        G[x][y] += sz * fr.dot(pa, pb)
                        Q[x][y] += sz * (2 * fr.dot(pa, pb) - fr.dot(pa, fr.nop(s, pb)))
                    if pa and nb in Cv:
                        Q[x][y] += sz * fr.dot(pa, Cv[nb])
                    if pb and na in Cv:
                        Q[x][y] += sz * fr.dot(pb, Cv[na])
    # pair terms
    for (sa, sb, sz) in pair_classes(be, mode):
        fr = Frame2(be, frame_labels(sa, sb))
        for la in letters:
            for lb in letters:
                a, b = (sa, la), (sb, lb)
                compa = {nm: fr.psi(nm, a) for nm in names}
                compb = {nm: fr.psi(nm, b) for nm in names}
                Ec = {nm: fr.Einv(fr.T(a, compb[nm])) for nm in names if compb[nm]}
                for x, na in enumerate(names):
                    pa = compa[na]
                    if not pa:
                        continue
                    Tb_pa = fr.T(b, pa)
                    smb_pa = fr.sm(sb, lb, pa)
                    for y, nb in enumerate(names):
                        pb = compb[nb]
                        if not pb:
                            continue
                        coef = be.R(1, 2) - (1 if sa == sb else 0)
                        Q[x][y] += sz * coef * fr.dot(smb_pa, fr.sm(sa, la, pb))
                        Q[x][y] -= sz * fr.dot(Tb_pa, Ec[nb])
    # G0 terms
    if mode == 'sym':
        fr = Frame2(be, [])
        ph = fr.phi('G0')
        g = names.index('G0')
        if ph:
            G[g][g] += fr.dot(ph, ph)
            Q[g][g] += 2 * fr.dot(ph, ph)
            Uph = fr.U(ph)
            if Uph:
                Q[g][g] -= fr.dot(Uph, fr.Einv(Uph))
    keep = [x for x in range(n) if G[x][x] != 0]
    return [names[x] for x in keep], [[G[x][y] for y in keep] for x in keep], [[Q[x][y] for y in keep] for x in keep]
