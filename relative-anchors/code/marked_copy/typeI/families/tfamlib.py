"""The orbit calculus of tlib.py (type I blocks with a letter 1, one marked copy) with the sizes symbolic: the Dicke
count q0 = j - tau - i, t = m - q0 (so m = q0 + t unpaired copies of U'), the irrep index i and tau' (the unmarked
copies of T), each fixed or symbolic; sigma is always symbolic.  Entries are rational functions in flintfield.Field
(python-flint fmpq_mpoly).  For both cases k = tau' + q0 + t + 2i + 1; case MT: tau = tau' + 1, j = q0 + tau' + 1 + i;
case MU: tau = tau', j = q0 + tau' + i.

Formal continuation (as famlib.py is for mlib.py):
 - normalisation: the numeric form of tlib.py at a member equals N0 times this one, N0 = 2^i C(m, q0 + sh) > 0
   (norm0; sh = -q0 if q0 is fixed and negative, sh = t if t is fixed and negative, else 0);
 - Dicke weights 2^{-pD} C(m - r, q0 + d) / C(m, q0 + sh) and the raising and lowering factors (q0 + d + 1),
   (t - r - d + 1) are polynomial or rational in the sizes; range checks are kept only for fixed sizes.  A key out of
   the range of a symbolic size (a phantom) has weight 0 at that member and the factors that would take it back into
   the range vanish there, so the part inside the range evolves as under the true operators;
 - the energy sector of a key is (kind, bT, delta): kind 'N' if copy 0 carries a 2, 'B' if it is Boolean; bT the
   Boolean T'-copies of the frame; delta such that the letters 1 (unmarked for 'N', all for 'B') number q0 + i + delta.
   Then n = m + 2i + bT, kap2 = tau' - bT, and every quantity of Lemma E' at the node i' = i + e is a polynomial in
   (q0, t, i, tau', sigma);
 - Newton nodes: e in -2..2, each condition of genuineness decided exactly when it involves only fixed sizes and
   taken as true when it involves a symbolic size (true for large members).  At a member where fewer nodes are
   genuine the extra nodes do not change the result on the part inside the range (it has no component there; at a
   node genuine on one side only the two-sided formulas reduce to the one-sided ones since mu = 0 there);
 - every denominator introduced on the way (node values, differences of node values, weight factors) is recorded in
   C.dens; the certificates check that each is one-signed after the shifts, hence nonzero at every member and sigma > 0.
Hence at every member and every sigma > 0 the family form, evaluated there, is the true reduced form divided by N0
(checked against tlib.py at sample points by the certificates)."""
import itertools
from fractions import Fraction as Fr
from math import comb
import flint
from flintfield import Field
import tlib

SIZES = ('q0', 't', 'i', 'tp')
SMIN = {'q0': 4, 't': 4, 'i': 2, 'tp': 2}
NAMES = tlib.NAMES
LETTERS = (1, 2)


def KF(K, x):
    x = Fr(x)
    return K(flint.fmpq(x.numerator, x.denominator))


class NonUniform(Exception):
    pass


def member(case, q0, t, i, tp):
    """(k, j, tau) of the member with these sizes"""
    k = tp + q0 + t + 2 * i + 1
    tau = tp + 1 if case == 'MT' else tp
    return k, q0 + tau + i, tau


class FamCalc:
    def __init__(self, case, mode, q0=None, t=None, i=None, tp=None, bscale=1, smin=None):
        assert case in ('MT', 'MU') and mode in ('sym', 'anti')
        self.case, self.mode = case, mode
        self.fixed = {'q0': q0, 't': t, 'i': i, 'tp': tp}
        self.smin = dict(SMIN)
        if smin:
            self.smin.update(smin)
        self.syms = [nm for nm in SIZES if self.fixed[nm] is None] + ['s']
        self.K = Field(self.syms)
        g = dict(zip(self.syms, self.K.gens))
        val = lambda nm: g[nm] if self.fixed[nm] is None else self.K(self.fixed[nm])
        self.q0, self.t, self.i, self.tp = (val(nm) for nm in SIZES)
        self.sig = g['s']
        self.one, self.zero = self.K(1), self.K(0)
        self.m = self.q0 + self.t
        self.up = self.m + 2 * self.i
        self.bscale = KF(self.K, bscale)
        sh = 0
        if q0 is not None and q0 < 0:
            sh = -q0
        if t is not None and t < 0:
            sh = t
        self.sh = sh
        self._einv = {}
        self.dens = {}
        self.stats = {'sectors': set(), 'dmin': 0, 'dmax': 0}

    # ---- sizes and uniform conditions
    def mn(self, nm):
        if nm == 'm':
            return self.mn('q0') + self.mn('t')
        return self.smin[nm] if self.fixed[nm] is None else self.fixed[nm]

    def symbolic(self, nm):
        if nm == 'm':
            return self.fixed['q0'] is None or self.fixed['t'] is None
        return self.fixed[nm] is None

    def ge(self, nm, c):
        """size nm >= c at every member (raises NonUniform if it holds only at some)"""
        if not self.symbolic(nm):
            return self.mn(nm) >= c
        if self.mn(nm) >= c:
            return True
        raise NonUniform(f'{nm} >= {c} holds only for some members')

    def holds(self, coeffs, c):
        """sum coeffs[nm] * nm >= c for large members: True if a symbolic size has a positive coefficient, otherwise
        decided exactly (the coefficients are nonnegative)"""
        if any(self.symbolic(nm) for nm, a in coeffs.items() if a > 0):
            return True
        return sum(a * self.fixed[nm] for nm, a in coeffs.items()) >= c

    def need(self, x, why):
        """record a denominator (a field element whose numerator must be one-signed after the shifts)"""
        if x == 0:
            raise ZeroDivisionError(why)
        key = (str(x.n), str(x.d))
        if key not in self.dens:
            self.dens[key] = (x, why)
        return x

    def member(self, vals):
        return member(self.case, vals['q0'], vals['t'], vals['i'], vals['tp'])

    def minpoint(self):
        return {nm: self.mn(nm) for nm in SIZES}

    # ---- E^-inverse data per sector (kind, bT, delta)
    def einv_data(self, sig):
        if sig in self._einv:
            return self._einv[sig]
        kind, bT, dl = sig
        self.stats['sectors'].add(sig)
        K, s = self.K, self.sig
        q0, t, i, tp = self.q0, self.t, self.i, self.tp
        half = KF(K, Fr(1, 2))
        K2 = tp - bT
        win = range(-2, 3)
        # n - c and c (the letters 1 counted as in Lemma E'), as polynomials
        nmc = t + i + (bT - dl)
        cc = q0 + i + dl

        def ok_i(e):
            return self.holds({'i': 1}, -e)

        if kind == 'N':
            nodes = [e for e in win if ok_i(e) and self.holds({'q0': 1}, e - dl) and self.holds({'t': 1}, e - bT + dl)]
            lam = lambda e: (t + (bT - dl - e)) * (q0 + (dl + 1 - e))
            vals = {}
            for e in nodes:
                vals[e] = self.need(K2 + s + KF(K, Fr(3, 2)) * nmc + half * cc + half * lam(e), ('N value', sig, e))
            data = {'N': self._newton(lambda e: self.one / vals[e], nodes, [lam(e) for e in nodes], sig)}
        else:
            inA = lambda e: ok_i(e) and self.holds({'q0': 1}, e - dl) and self.holds({'t': 1}, e - bT + dl)
            inB = lambda e: ok_i(e) and self.holds({'q0': 1}, e - dl + 1) and self.holds({'t': 1}, e - bT + dl - 1)
            nA = [e for e in win if inA(e)]
            nB = [e for e in win if inB(e)]
            lamA = lambda e: (t + (bT - dl - e)) * (q0 + (dl + 1 - e))
            lamB = lambda e: (t + (bT - dl + 1 - e)) * (q0 + (dl - e))
            vals = {}
            for e in sorted(set(nA) | set(nB)):
                mu = lamB(e)
                nu = (t + (bT - dl + 2 - e)) * (q0 + (dl - 1 - e))
                a = K2 + 2 * s + 2 * nmc + half * mu
                ee = K2 + half * s + 2 * (nmc + 1) + half * nu
                if e in nA and e in nB:
                    D = a * ee - s * mu * KF(K, Fr(1, 4))
                    vals[e] = ('both', a, ee, self.need(D, ('B Delta', sig, e)))
                elif e in nA:
                    vals[e] = ('A', self.need(a, ('B a', sig, e)), ee, a)
                else:
                    vals[e] = ('B', a, self.need(ee, ('B e', sig, e)), ee)

            def f1(e):
                kd, a, ee, D = vals[e]
                return ee / D if kd == 'both' else self.one / a

            def f2(e):
                kd, a, ee, D = vals[e]
                return (-half * s) / D if kd == 'both' else self.zero

            def f3(e):
                kd, a, ee, D = vals[e]
                return (s * a) / D if kd == 'both' else s / ee
            data = {'A1': self._newton(f1, nA, [lamA(e) for e in nA], sig),
                    'B2': self._newton(f2, nB, [lamB(e) for e in nB], sig),
                    'B3': self._newton(f3, nB, [lamB(e) for e in nB], sig)}
        self._einv[sig] = data
        return data

    def _newton(self, fval, nodes, lamv, sig):
        coefs, diag = [], []
        for n_, e in enumerate(nodes):
            new = [fval(e)]
            for l in range(1, n_ + 1):
                dd = self.need(lamv[n_] - lamv[n_ - l], ('node difference', sig, nodes[n_], nodes[n_ - l]))
                new.append((new[l - 1] - diag[l - 1]) / dd)
            diag = new
            coefs.append(new[n_])
        return lamv, coefs

    # ---- Dicke weights (normalised by C(m, q0 + sh) 2^i)
    def weight(self, pD, r, d):
        K = self.K
        m = self.m
        q0 = self.q0 + self.sh
        d = d - self.sh
        w = KF(K, Fr(1, 2 ** pD))
        if d >= 0:
            for s_ in range(1, d + 1):
                w = w / self.need(q0 + s_, 'weight q0 + s')
        else:
            for s_ in range(0, -d):
                w = w * (q0 - s_)
        s2 = r + d
        if s2 >= 0:
            for s_ in range(0, s2):
                w = w * (m - q0 - s_)
        else:
            for s_ in range(1, -s2 + 1):
                w = w / self.need(m - q0 + s_, 'weight t + s')
        for s_ in range(0, r):
            w = w / self.need(m - s_, 'weight m - s')
        return w

    def inrange(self, d, r):
        """is the Dicke offset d (r distinguished copies of R) in range?  exact for fixed sizes, True (phantoms
        allowed) for symbolic ones"""
        if self.fixed['q0'] is not None and self.fixed['q0'] + d < 0:
            return False
        if self.fixed['t'] is not None and d > self.fixed['t'] - r:
            return False
        return True

    def valid(self):
        if self.mode == 'anti' and not self.ge('tp', 2):
            return False
        for nm in SIZES:
            if not self.symbolic(nm) and self.fixed[nm] < 0 and nm in ('i', 'tp'):
                return False
        if not self.symbolic('m') and self.mn('m') < 0:
            return False
        if self.case == 'MU' and not self.ge('tp', 1):
            return False
        return True


def frame_labels(*sites):
    return tlib.frame_labels(*sites)


class Frame:
    def __init__(self, C, labels):
        self.C = C
        self.labels = labels
        assert labels[0] == 'M'
        self.pos = {l: n for n, l in enumerate(labels)}
        self.pairs = sorted({l[1] for l in labels if l[0] == 'P'})
        self.rs = [l for l in labels if l[0] == 'R']
        self.ts = [l for l in labels if l[0] == 'T']
        self.pD = len(self.pairs)
        self.r = len(self.rs)
        self.upos = [self.pos[l] for l in labels if l[0] in 'PR']
        self.tpos = [self.pos[l] for l in self.ts]
        self.unm = self.upos + self.tpos
        K = C.K
        self.h, self.th = KF(K, Fr(1, 2)), KF(K, Fr(3, 2))
        self.mh = KF(K, Fr(-1, 2))

    def inrange(self, d):
        return self.C.inrange(d, self.r)

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

    def note(self, d):
        st = self.C.stats
        st['dmin'] = min(st['dmin'], d)
        st['dmax'] = max(st['dmax'], d)

    def jUp(self, v):
        out = {}
        q0 = self.C.q0
        for (L, d), c in v.items():
            for s in self.unm:
                if L[s] == 0:
                    self.addto(out, (L[:s] + (1,) + L[s + 1:], d), c)
            if self.inrange(d + 1):
                self.note(d + 1)
                self.addto(out, (L, d + 1), c * (q0 + (d + 1)))
        return out

    def jUm(self, v):
        out = {}
        t = self.C.t
        for (L, d), c in v.items():
            for s in self.unm:
                if L[s] == 1:
                    self.addto(out, (L[:s] + (0,) + L[s + 1:], d), c)
            if self.inrange(d - 1):
                self.note(d - 1)
                self.addto(out, (L, d - 1), c * (t - (self.r + d - 1)))
        return out

    def flipM(self, v, frm, to):
        return {(((to,) + L[1:]), d): c for (L, d), c in v.items() if L[0] == frm}

    def jp(self, v):
        return self.add(self.jUp(v), self.flipM(v, 0, 1))

    def jm(self, v):
        return self.add(self.jUm(v), self.flipM(v, 1, 0))

    def sm(self, lab, l, v):
        s = self.pos[lab]
        return {(L[:s] + (0,) + L[s + 1:], d): c for (L, d), c in v.items() if L[s] == l}

    def nop(self, lab, v):
        s = self.pos[lab]
        return {key: c for key, c in v.items() if key[0][s] != 0}

    def Nop(self, lab, v):
        s = self.pos[lab]
        f = {0: self.th, 1: self.h, 2: self.C.one}
        return {(L, d): c * f[L[s]] for (L, d), c in v.items()}

    def t(self, a, v):
        lab, l = a
        out = self.scale(self.sm(lab, l, self.jp(v)), self.h)
        if l == 1:
            return self.add(out, self.Nop(lab, v))
        s = self.pos[lab]
        for (L, d), c in v.items():
            if L[s] == 2:
                self.addto(out, (L[:s] + (1,) + L[s + 1:], d), c * self.mh)
        return out

    def Pi(self, v):
        s = self.C.sig
        return {(L, d): (c * s if L[0] == 1 else c) for (L, d), c in v.items()}

    def DU(self, v):
        return self.jUm(self.jUp(v))

    def sector(self, key):
        L, d = key
        bT = sum(1 for s in self.tpos if L[s] != 2)
        f1 = sum(1 for s in self.unm if L[s] == 1)
        dl = f1 - self.pD + d
        if L[0] == 2:
            return ('N', bT, dl)
        return ('B', bT, dl + L[0])

    def einv(self, v):
        groups = {}
        for key, x in v.items():
            groups.setdefault(self.sector(key), {})[key] = x
        out = {}
        for sk, w in groups.items():
            data = self.C.einv_data(sk)
            if sk[0] == 'N':
                self.add(out, self.newton(w, *data['N']))
            else:
                vA = {key: x for key, x in w.items() if key[0][0] == 0}
                vB = {key: x for key, x in w.items() if key[0][0] == 1}
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

    def dot(self, a, b):
        if len(a) > len(b):
            a, b = b, a
        s = self.C.zero
        for key, x in a.items():
            y = b.get(key)
            if y is not None:
                s = s + x * y * self.C.weight(self.pD, self.r, key[1])
        return s

    def _expand(self, fixed, qtot, sign=1):
        C = self.C
        free_pairs = [p for p in self.pairs if ('P%sa' % p) not in fixed]
        free_rs = [r for r in self.rs if r not in fixed]
        if 'M' not in fixed:
            assert C.case == 'MT'
        out = {}
        for sing in itertools.product((0, 1), repeat=len(free_pairs)):
            for rl in itertools.product((0, 1), repeat=len(free_rs)):
                L = [None] * len(self.labels)
                coef = sign
                L[0] = 2
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
                d = qtot - sum(rl)
                if not self.inrange(d):
                    continue
                self.note(d)
                self.addto(out, (tuple(L), d), self.C.one * coef)
        return out

    # the twisted vectors: exactly those of tlib.Frame (same code, field coefficients)
    psi = tlib.Frame.psi
    phi = tlib.Frame.phi


def orbits(C):
    """H-orbits of row sites and of ordered pairs of row sites (label, size), as tlib.orbits, sizes symbolic"""
    m, i, tp = C.m, C.i, C.tp
    if C.mode == 'anti':
        two = C.K(2)
        return [('Ta', two)], [('Ta', 'Ta', two), ('Ta', 'Tb', two)]
    sr = [('M', C.one)]
    if C.ge('m', 1):
        sr.append(('R0', m))
    if C.ge('i', 1):
        sr.append(('P0a', 2 * i))
    if C.ge('tp', 1):
        sr.append(('Ta', tp))
    pr = [('M', 'M', C.one)]
    if C.ge('m', 1):
        pr += [('M', 'R0', m), ('R0', 'M', m), ('R0', 'R0', m)]
    if C.ge('m', 2):
        pr.append(('R0', 'R1', m * (m - 1)))
    if C.ge('i', 1):
        pr += [('M', 'P0a', 2 * i), ('P0a', 'M', 2 * i), ('P0a', 'P0a', 2 * i), ('P0a', 'P0b', 2 * i)]
        if C.ge('m', 1):
            pr += [('R0', 'P0a', 2 * i * m), ('P0a', 'R0', 2 * i * m)]
    if C.ge('i', 2):
        pr.append(('P0a', 'P1a', 2 * i * (2 * i - 2)))
    if C.ge('tp', 1):
        pr += [('M', 'Ta', tp), ('Ta', 'M', tp), ('Ta', 'Ta', tp)]
        if C.ge('m', 1):
            pr += [('R0', 'Ta', m * tp), ('Ta', 'R0', m * tp)]
        if C.ge('i', 1):
            pr += [('P0a', 'Ta', 2 * i * tp), ('Ta', 'P0a', 2 * i * tp)]
    if C.ge('tp', 2):
        pr.append(('Ta', 'Tb', tp * (tp - 1)))
    return sr, pr


def present_at(case, mode, vals):
    """names with a nonzero twisted vector at the member with these sizes (tlib.norm_of)"""
    k, j, tau = member(case, vals['q0'], vals['t'], vals['i'], vals['tp'])
    C = tlib.Calc(k, j, tau, vals['i'], case, mode, sigma=1)
    if not C.valid():
        return []
    return [nm for nm in NAMES[(case, mode)] if tlib.norm_of(C, nm) != 0]


def present_names(C):
    """the twisted vectors present on the family: decided at the minimal member; every condition of tlib.norm_of is a
    lower bound on a nonnegative combination of the sizes, so presence is monotone, and it must agree with the member
    where every symbolic size is increased by 100 (otherwise the family is not uniform)"""
    lo = C.minpoint()
    hi = {nm: lo[nm] + (100 if C.symbolic(nm) else 0) for nm in SIZES}
    a, b = present_at(C.case, C.mode, lo), present_at(C.case, C.mode, hi)
    if a != b:
        raise NonUniform(f'presence of vectors differs: {a} at the minimum, {b} far out')
    return a


def reduced(C):
    """(names, Q) over the field: the reduced form on the present twisted vectors, divided by N0"""
    if not C.valid():
        return [], []
    names = present_names(C)
    if not names:
        return [], []
    sr, pr = orbits(C)
    sig = C.sig
    n = len(names)
    ix = {nm: x for x, nm in enumerate(names)}
    Q = [[C.zero] * n for _ in range(n)]
    b = C.bscale
    s_of = lambda a: sig if (a[0] == 'M' and a[1] == 1) else C.one
    for (st, sz) in sr:
        fr = Frame(C, frame_labels(st))
        phis = {nm: fr.phi(nm) for nm in names}
        for l in LETTERS:
            a = (st, l)
            comps = {nm: fr.psi(nm, a) for nm in names}
            if not any(comps.values()):
                continue
            sa = s_of(a)
            Cv = {}
            for nm in names:
                ph = phis[nm]
                if ph:
                    v = fr.scale(fr.Pi(fr.sm(st, l, ph)), -sa)
                    fr.add(v, fr.scale(fr.jm(fr.einv(fr.t(a, ph))), sa * b))
                    Cv[nm] = v
            for na in names:
                pa = comps[na]
                if not pa:
                    continue
                x = ix[na]
                Pipa = fr.Pi(pa)
                for nb in names:
                    y = ix[nb]
                    pb = comps[nb]
                    if pb:
                        Q[x][y] += sz * sa * (2 * fr.dot(Pipa, pb) - fr.dot(Pipa, fr.nop(st, pb)))
                    if nb in Cv:
                        v0 = sz * fr.dot(pa, Cv[nb])
                        Q[x][y] += v0
                        Q[y][x] += v0
    half = KF(C.K, Fr(1, 2))
    for (sa_, sb_, sz) in pr:
        fr = Frame(C, frame_labels(sa_, sb_))
        for la in LETTERS:
            for lb in LETTERS:
                a, bb = (sa_, la), (sb_, lb)
                compa = {nm: fr.psi(nm, a) for nm in names}
                compb = {nm: fr.psi(nm, bb) for nm in names}
                if not any(compa.values()) or not any(compb.values()):
                    continue
                sab = s_of(a) * s_of(bb)
                coef = half - (1 if sa_ == sb_ else 0)
                Ecache, smA = {}, {}
                for nb in names:
                    pb = compb[nb]
                    if pb:
                        Ecache[nb] = fr.einv(fr.t(a, pb))
                        smA[nb] = fr.Pi(fr.sm(sa_, la, pb))
                for na in names:
                    pa = compa[na]
                    if not pa:
                        continue
                    x = ix[na]
                    Tb_pa = fr.t(bb, pa)
                    smb_pa = fr.sm(sb_, lb, pa)
                    for nb in Ecache:
                        y = ix[nb]
                        Q[x][y] += sz * (coef * sab * fr.dot(smb_pa, smA[nb]) - b * sab * fr.dot(Tb_pa, Ecache[nb]))
    if C.mode == 'sym':
        fr = Frame(C, ['M'])
        ph = {nm: fr.phi(nm) for nm in names}
        for na in names:
            for nb in names:
                if ph[na] and ph[nb]:
                    x, y = ix[na], ix[nb]
                    Q[x][y] += 2 * fr.dot(fr.Pi(ph[na]), ph[nb])
                    Ua, Ub = fr.jp(ph[na]), fr.jp(ph[nb])
                    if Ua and Ub:
                        Q[x][y] += -b * fr.dot(Ua, fr.einv(Ub))
    return names, Q


def norm0(C, vals):
    """N0 = 2^i C(m, q0 + sh) at the member with these sizes: tlib's numeric form there equals N0 times this one"""
    m = vals['q0'] + vals['t']
    return (2 ** vals['i']) * comb(m, vals['q0'] + C.sh)


def kernel_vectors(C, names):
    """the twisted parts of p (x) H_j (tlib.kernel_vectors), with field coefficients"""
    return [(lab, {nm: C.K(cf) for nm, cf in v.items()}) for lab, v in tlib.kernel_vectors(C, names)]
