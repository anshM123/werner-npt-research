"""Orbit calculus for the type I torus blocks (with a letter 1) of (R_inf)_j with one marked copy.

Squared weights: copy 0 (label 'M') has sigma > 0, every other copy has 1.  T = the copies that carry one letter 2
(row or column), |T| = tau, 1 <= tau <= j; U = the other copies.  Two cases:
  case 'MT': M in T.  T' = T - {M} (tau' = tau - 1 unmarked T-copies), U' = U (u' = k - tau copies);
  case 'MU': M in U.  T' = T (tau' = tau),                            U' = U - {M} (u' = k - tau - 1 copies).
Symmetry S_{T'} x S_{U'}.  T'-modes: 'sym' (S_{T'}-invariant) and 'anti' (sign on the pair (Ta, Tb) of T', the
irrep (tau' - 1, 1)).  U'-twist: H_i = (S_2 wr S_i) x S_m on U' (i pairs, m = u' - 2i rest copies R), character
chi_i = sign on each pair swap.  q0 = j - tau - i is the base Dicke count on R.

Rescaled frame (each state multiplied by mu = prod_t w_t^{N_t}, N_t the letters 1 at copy t in row and column), units
2^Z (Z = k - j); with Pi = diag(sigma^{[letter 1 at M]}), s_a^{l1} = sigma^{[row a = 1@M]}, jb^{+-} unweighted:
  2^Z Q(Psi) = sum_{a,b} (1/2 - [r_a = r_b]) s_a^{l1} s_b^{l1} <sig_b^- psi_a, Pi sig_a^- psi_b>
             + sum_a s_a^{l1} <psi_a, Pi (2 - n_{r_a}) psi_a>
             - sum_{a,b} s_a^{l1} s_b^{l1} <t_b psi_a, E^^{-1} t_a psi_b>
             - 2 sum_a s_a^{l1} ( <psi_a, Pi sig_a^- phi> - <jb^+ psi_a, E^^{-1} t_a phi> )
             + 2 <phi, Pi phi> - <jb^+ phi, E^^{-1} jb^+ phi>,
  t_(r,1) = 1/2 sig_(r,1)^- jb^+ + N_r (N = 3/2, 1/2, 1 on the letters 0, 1, 2),
  t_(t,2) = 1/2 sig_(t,2)^- jb^+ - 1/2 |1><2|_t,
  E^ = D_w^{-1} E D_w^{-1}, E = sum_t s_t N_t + 1/2 J^- J^+.
On strings with a fixed pattern of letters 2: n = number of unmarked copies without a letter 2, kap2 = number of
unmarked letters 2 (+ sigma if M carries a 2), c' = letters 1 on unmarked copies.
  M carries 2:   E^ = E = kap2 + 3/2 (n - c') + 1/2 c' + 1/2 DU, DU = jb_U^- jb_U^+ (Johnson operator on n copies).
  M Boolean, c ones in all:  A = letter 0 at M, B = letter 1 at M;
     E^ = [[kap2 + 2 sigma + 2(n - c) + 1/2 J+J-,  1/2 J+ sig_M^-], [1/2 sig_M^+ J-,  (kap2 + 2(n - c + 1) + sigma/2
            + 1/2 J+J-)/sigma]]
  inverted on each S_n-irrep by the 2 x 2 formula (Lemma E of marked_copy/NOTES.md with kap2 added); E^-inverse is
  applied by Newton interpolation in DU over the Johnson nodes of a window, and the run asserts that the remainder
  vanishes (so the window is exact on the vector).
Frame vectors: dict (L, d) -> coefficient; L = letters on the frame labels ('M' first, then U'-labels 'P0a','P0b',
'P1a','P1b','R0','R1', then T'-labels 'Ta','Tb'); undistinguished pairs carry singlets, undistinguished T'-copies the
letter 2, the rest R - frame carries Dicke(q0 + d).
Modes: Calc(..., sigma=Fraction) exact numeric (fmpq); Calc(..., sigma=None) symbolic in sigma (fmpq_poly
coefficients; every E^-inverse returns pieces (key, vector) meaning vector / D[key]).
bscale multiplies every term with an inverse E^-1 (the B-part), as in j1lib.py.
"""
import itertools
from fractions import Fraction as Fr
from math import comb
import flint
from ufield import RF

fq = flint.fmpq
PL = flint.fmpq_poly


def tofq(x):
    if isinstance(x, Fr):
        return fq(x.numerator, x.denominator)
    return fq(x)


NAMES = {
    ('MT', 'sym'): ['A0', 'A1', 'B00', 'B01', 'B10', 'B11', 'C', 'D0', 'D1', 'CM', 'DM0', 'DM1', 'G'],
    ('MT', 'anti'): ['Ca', 'Da0', 'Da1'],
    ('MU', 'sym'): ['O0', 'O1', 'A00', 'A01', 'A10', 'A11', 'B000', 'B001', 'B010', 'B011', 'B100', 'B101', 'B110',
                    'B111', 'C0', 'C1', 'D00', 'D01', 'D10', 'D11', 'G0', 'G1'],
    ('MU', 'anti'): ['Ca0', 'Ca1', 'Da00', 'Da01', 'Da10', 'Da11'],
}


class Calc:
    def __init__(self, k, j, tau, i, case, mode, sigma=None, bscale=1, window=2):
        assert case in ('MT', 'MU') and mode in ('sym', 'anti')
        self.k, self.j, self.tau, self.i, self.case, self.mode = k, j, tau, i, case, mode
        self.taup = tau - 1 if case == 'MT' else tau
        self.up = k - tau if case == 'MT' else k - tau - 1
        self.m = self.up - 2 * i
        self.q0 = j - tau - i
        self.sym = sigma is None
        if self.sym:
            self.sig = PL([0, 1])
            self.one = PL([1])
        else:
            self.sig = tofq(sigma)
            self.one = fq(1)
        self.zero = self.one * 0
        self.bscale = tofq(bscale)
        self.window = window
        self._einv = {}
        self.D = {}

    def valid(self):
        if self.taup < 0 or self.up < 0 or self.m < 0 or self.i < 0:
            return False
        if self.mode == 'anti' and self.taup < 2:
            return False
        return True

    def winnodes(self, nodes):
        return [x for x in nodes if self.i - self.window <= x <= self.i + self.window]

    # ---------------- E^-inverse data per sector signature ----------------
    def einv_data(self, sigkey):
        """sigkey = ('N', n, cp, kap2) for M carrying a 2 (kap2 = unmarked 2s; sigma added), or ('B', n, c, kap2) for
        M Boolean with c letters 1 in all.  Returns Newton data (lambdas, coefficients) per kind, and D."""
        if sigkey in self._einv:
            return self._einv[sigkey]
        s = self.sig
        kind, n, c, kap2 = sigkey
        K2 = s * 0 + kap2
        data = {}
        if kind == 'N':
            cp = c
            IA = list(range(0, min(cp, n - cp) + 1)) if 0 <= cp <= n else []
            nA = self.winnodes(IA)
            vals = {}
            for ip in nA:
                lam = (n - cp - ip) * (cp + 1 - ip)
                vals[ip] = K2 + s + tofq(Fr(3, 2) * (n - cp) + Fr(cp, 2) + Fr(lam, 2))
            Dc = self.one
            if self.sym:
                for ip in nA:
                    Dc = Dc * vals[ip]

            def fv(ip):
                return (self.one, vals[ip])
            data['N'] = self._newton(fv, nA, [(n - cp - ip) * (cp + 1 - ip) for ip in nA], Dc)
        else:
            IA = list(range(0, min(c, n - c) + 1)) if 0 <= c <= n else []
            IB = list(range(0, min(c - 1, n - c + 1) + 1)) if 0 <= c - 1 <= n else []
            vals = {}
            for ip in sorted(set(IA) | set(IB)):
                mu = (n - c + 1 - ip) * (c - ip)
                nu = (n - c + 2 - ip) * (c - 1 - ip)
                a = K2 + 2 * s + tofq(Fr(2 * (n - c)) + Fr(mu, 2))
                e = K2 + s * fq(1, 2) + tofq(Fr(2 * (n - c + 1)) + Fr(nu, 2))
                if ip in IA and ip in IB:
                    vals[ip] = ('both', a, e, a * e - s * fq(mu, 4))
                elif ip in IA:
                    vals[ip] = ('A', a, e, a)
                else:
                    vals[ip] = ('B', a, e, e)
            nA = self.winnodes(IA)
            nB = self.winnodes(IB)
            used = sorted(set(nA) | set(nB))
            Dc = self.one
            if self.sym:
                for ip in used:
                    Dc = Dc * vals[ip][3]

            def f1(ip):
                t, a, e, Dd = vals[ip]
                return (e, Dd) if t == 'both' else (self.one, a)

            def f2(ip):
                t, a, e, Dd = vals[ip]
                return ((-s * fq(1, 2)), Dd) if t == 'both' else (self.zero, self.one)

            def f3(ip):
                t, a, e, Dd = vals[ip]
                return (s * a, Dd) if t == 'both' else (s, e)
            lamA = [(n - c - ip) * (c + 1 - ip) for ip in nA]
            lamB = [(n - c + 1 - ip) * (c - ip) for ip in nB]
            data['A1'] = self._newton(f1, nA, lamA, Dc)
            data['B2'] = self._newton(f2, nB, lamB, Dc)
            data['B3'] = self._newton(f3, nB, lamB, Dc)
        data['D'] = Dc
        self._einv[sigkey] = data
        self.D[sigkey] = Dc
        return data

    def _newton(self, fval, nodes, lamv, Dc):
        coefs = []
        diag = []
        for n_, ip in enumerate(nodes):
            num, den = fval(ip)
            fv = RF(num, den) if self.sym else num / den
            new = [fv]
            for l in range(1, n_ + 1):
                new.append((new[l - 1] - diag[l - 1]) / (lamv[n_] - lamv[n_ - l]))
            diag = new
            coefs.append(new[n_])
        if self.sym:
            out = []
            for cf in coefs:
                x = cf * RF(Dc)
                assert x.d == PL([1]), 'Newton coefficient times D is not a polynomial'
                out.append(x.n)
            coefs = out
        return lamv, coefs


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
    # order: M, U'-labels, T'-labels
    return ['M'] + [l for l in labs[1:] if l[0] in 'PR'] + [l for l in labs[1:] if l[0] == 'T']


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
        self.nprime = C.m - len(self.rs)
        self.undist = C.i - self.pD
        self.upos = [self.pos[l] for l in labels if l[0] in 'PR']      # U'-positions (always Boolean)
        self.tpos = [self.pos[l] for l in self.ts]                      # T'-positions (Boolean unless letter 2)
        self.unm = self.upos + self.tpos                               # all unmarked distinguished positions

    def inrange(self, d):
        q = self.C.q0 + d
        return 0 <= q <= self.nprime

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

    # ---- unweighted raising / lowering on the unmarked Boolean copies (letter 2 inert; singlets killed)
    def jUp(self, v):
        out = {}
        q0 = self.C.q0
        for (L, d), c in v.items():
            for s in self.unm:
                if L[s] == 0:
                    self.addto(out, (L[:s] + (1,) + L[s + 1:], d), c)
            if self.inrange(d + 1):
                self.addto(out, (L, d + 1), c * (q0 + d + 1))
        return out

    def jUm(self, v):
        out = {}
        q0 = self.C.q0
        for (L, d), c in v.items():
            for s in self.unm:
                if L[s] == 1:
                    self.addto(out, (L[:s] + (0,) + L[s + 1:], d), c)
            if self.inrange(d - 1):
                self.addto(out, (L, d - 1), c * (self.nprime - (q0 + d) + 1))
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
        f = {0: fq(3, 2), 1: fq(1, 2), 2: fq(1)}
        return {(L, d): c * f[L[s]] for (L, d), c in v.items()}

    def t(self, a, v):
        lab, l = a
        out = self.scale(self.sm(lab, l, self.jp(v)), fq(1, 2))
        if l == 1:
            return self.add(out, self.Nop(lab, v))
        s = self.pos[lab]
        for (L, d), c in v.items():
            if L[s] == 2:
                self.addto(out, (L[:s] + (1,) + L[s + 1:], d), c * fq(-1, 2))
        return out

    def Pi(self, v):
        s = self.C.sig
        return {(L, d): (c * s if L[0] == 1 else c) for (L, d), c in v.items()}

    def DU(self, v):
        return self.jUm(self.jUp(v))

    def sector(self, key):
        L, d = key
        C = self.C
        bT = sum(1 for s in self.tpos if L[s] != 2)
        n = C.up + bT
        kap2 = C.taup - bT
        cp = sum(1 for s in self.unm if L[s] == 1) + self.undist + C.q0 + d
        if L[0] == 2:
            return ('N', n, cp, kap2)
        return ('B', n, cp + L[0], kap2)

    def einv(self, v):
        """E^-inverse of v; returns a list of (sigkey, vector): in symbolic mode the true inverse is vector / D[sigkey]"""
        groups = {}
        for key, x in v.items():
            groups.setdefault(self.sector(key), {})[key] = x
        res = []
        for sk, w in groups.items():
            data = self.C.einv_data(sk)
            out = {}
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
            if out:
                res.append((sk, out))
        return res

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
        if w:
            raise AssertionError('Newton window not exact (remainder nonzero)')
        return out

    def dot(self, a, b):
        if len(a) > len(b):
            a, b = b, a
        s = self.C.zero
        q0, und = self.C.q0, self.undist
        for key, x in a.items():
            y = b.get(key)
            if y is not None:
                s = s + x * y * ((2 ** und) * comb(self.nprime, q0 + key[1]))
        return s

    # ---- twisted vectors
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
                self.addto(out, (tuple(L), d), self.C.one * coef)
        return out

    def psi(self, name, a):
        """component at the row a = (site label, letter) of the twisted vector name (frame vector)"""
        C = self.C
        s, l = a
        if C.case == 'MT':
            if C.mode == 'sym':
                if name in ('A0', 'A1'):
                    x = int(name[1])
                    return self._expand({s: x}, -x) if (s[0] == 'R' and l == 1) else {}
                if name[0] == 'B':
                    if s[0] != 'P' or l != 1:
                        return {}
                    x, y = int(name[1]), int(name[2])
                    first = s[2] == 'a'
                    partner = 'P%s%s' % (s[1], 'b' if first else 'a')
                    return self._expand({s: x, partner: y}, 1 - x - y, 1 if first else -1)
                if name == 'C':
                    return self._expand({}, 0) if (s[0] == 'T' and l == 1) else {}
                if name in ('D0', 'D1'):
                    z = int(name[1])
                    return self._expand({s: z}, 1 - z) if (s[0] == 'T' and l == 2) else {}
                if name == 'CM':
                    return self._expand({}, 0) if (s == 'M' and l == 1) else {}
                if name in ('DM0', 'DM1'):
                    z = int(name[2])
                    return self._expand({'M': z}, 1 - z) if (s == 'M' and l == 2) else {}
                return {}
            # anti
            if s not in ('Ta', 'Tb'):
                return {}
            sg = 1 if s == 'Ta' else -1
            if name == 'Ca':
                return self._expand({}, 0, sg) if l == 1 else {}
            if name in ('Da0', 'Da1'):
                z = int(name[2])
                return self._expand({s: z}, 1 - z, sg) if l == 2 else {}
            return {}
        # case MU
        if C.mode == 'sym':
            if name[0] == 'O':
                y0 = int(name[1])
                return self._expand({'M': y0}, -y0) if (s == 'M' and l == 1) else {}
            if name[0] == 'A':
                if s[0] != 'R' or l != 1:
                    return {}
                x, y0 = int(name[1]), int(name[2])
                return self._expand({s: x, 'M': y0}, -x - y0)
            if name[0] == 'B':
                if s[0] != 'P' or l != 1:
                    return {}
                x, y, y0 = int(name[1]), int(name[2]), int(name[3])
                first = s[2] == 'a'
                partner = 'P%s%s' % (s[1], 'b' if first else 'a')
                return self._expand({s: x, partner: y, 'M': y0}, 1 - x - y - y0, 1 if first else -1)
            if name[0] == 'C':
                y0 = int(name[1])
                return self._expand({'M': y0}, -y0) if (s[0] == 'T' and l == 1) else {}
            if name[0] == 'D':
                z, y0 = int(name[1]), int(name[2])
                return self._expand({s: z, 'M': y0}, 1 - z - y0) if (s[0] == 'T' and l == 2) else {}
            return {}
        if s not in ('Ta', 'Tb'):
            return {}
        sg = 1 if s == 'Ta' else -1
        if name[:2] == 'Ca':
            y0 = int(name[2])
            return self._expand({'M': y0}, -y0, sg) if l == 1 else {}
        if name[:2] == 'Da':
            z, y0 = int(name[2]), int(name[3])
            return self._expand({s: z, 'M': y0}, 1 - z - y0, sg) if l == 2 else {}
        return {}

    def phi(self, name):
        C = self.C
        if C.mode != 'sym' or name[0] != 'G':
            return {}
        if C.case == 'MT':
            return self._expand({}, 1)
        z0 = int(name[1])
        return self._expand({'M': z0}, 1 - z0)


def orbits(C):
    """H-orbits of row sites and of ordered pairs of row sites (label, size), for the mode of C"""
    m, i, tp = C.m, C.i, C.taup
    if C.mode == 'anti':
        return [('Ta', 2)], [('Ta', 'Ta', 2), ('Ta', 'Tb', 2)]
    sr = [('M', 1)]
    if m >= 1:
        sr.append(('R0', m))
    if i >= 1:
        sr.append(('P0a', 2 * i))
    if tp >= 1:
        sr.append(('Ta', tp))
    pr = [('M', 'M', 1)]
    if m >= 1:
        pr += [('M', 'R0', m), ('R0', 'M', m), ('R0', 'R0', m)]
    if m >= 2:
        pr.append(('R0', 'R1', m * (m - 1)))
    if i >= 1:
        pr += [('M', 'P0a', 2 * i), ('P0a', 'M', 2 * i), ('P0a', 'P0a', 2 * i), ('P0a', 'P0b', 2 * i)]
        if m >= 1:
            pr += [('R0', 'P0a', 2 * i * m), ('P0a', 'R0', 2 * i * m)]
    if i >= 2:
        pr.append(('P0a', 'P1a', 2 * i * (2 * i - 2)))
    if tp >= 1:
        pr += [('M', 'Ta', tp), ('Ta', 'M', tp), ('Ta', 'Ta', tp)]
        if m >= 1:
            pr += [('R0', 'Ta', m * tp), ('Ta', 'R0', m * tp)]
        if i >= 1:
            pr += [('P0a', 'Ta', 2 * i * tp), ('Ta', 'P0a', 2 * i * tp)]
    if tp >= 2:
        pr.append(('Ta', 'Tb', tp * (tp - 1)))
    return sr, pr


LETTERS = (1, 2)


def reduced(C):
    """(names, norms, Q) on the nonzero twisted vectors.  numeric mode: Q entries fmpq; symbolic mode: RF in sigma."""
    names_all = NAMES[(C.case, C.mode)]
    n = len(names_all)
    if not C.valid():
        return [], [], []
    sr, pr = orbits(C)
    j = C.j
    sig = C.sig
    Q0 = [[C.zero] * n for _ in range(n)]
    QE = {}
    G = [[fq(0)] * n for _ in range(n)]
    b = C.bscale

    def addE(x, y, pieces_dot):
        for sk, val in pieces_dot:
            if val == 0:
                continue
            M = QE.get(sk)
            if M is None:
                M = QE[sk] = [[C.zero] * n for _ in range(n)]
            M[x][y] += val

    def s_of(a):
        return sig if (a[0] == 'M' and a[1] == 1) else C.one

    # site terms
    for (st, sz) in sr:
        fr = Frame(C, frame_labels(st))
        phis = {nm: fr.phi(nm) for nm in names_all}
        for l in LETTERS:
            a = (st, l)
            comps = {nm: fr.psi(nm, a) for nm in names_all}
            if not any(comps.values()):
                continue
            sa = s_of(a)
            Cbase, CE = {}, {}
            for nm in names_all:
                ph = phis[nm]
                if ph:
                    Cbase[nm] = fr.scale(fr.Pi(fr.sm(st, l, ph)), -sa)
                    pieces = fr.einv(fr.t(a, ph))
                    CE[nm] = [(sk, fr.scale(fr.jm(vv), sa * b)) for sk, vv in pieces]
            for x, na in enumerate(names_all):
                pa = comps[na]
                if not pa:
                    continue
                Pipa = fr.Pi(pa)
                for y, nb in enumerate(names_all):
                    pb = comps[nb]
                    if pb:
                        G[x][y] += sz * fr.dot(pa, pb) if not C.sym else 0
                        Q0[x][y] += sz * sa * (2 * fr.dot(Pipa, pb) - fr.dot(Pipa, fr.nop(st, pb)))
                    if nb in Cbase:
                        v0 = sz * fr.dot(pa, Cbase[nb])
                        Q0[x][y] += v0
                        Q0[y][x] += v0
                        dd = [(sk, sz * fr.dot(pa, vv)) for sk, vv in CE[nb]]
                        addE(x, y, dd)
                        addE(y, x, dd)
    # pair terms
    for (sa_, sb_, sz) in pr:
        fr = Frame(C, frame_labels(sa_, sb_))
        for la in LETTERS:
            for lb in LETTERS:
                a, bb = (sa_, la), (sb_, lb)
                compa = {nm: fr.psi(nm, a) for nm in names_all}
                compb = {nm: fr.psi(nm, bb) for nm in names_all}
                if not any(compa.values()) or not any(compb.values()):
                    continue
                sab = s_of(a) * s_of(bb)
                coef = fq(1, 2) - (1 if sa_ == sb_ else 0)
                Ecache, smA = {}, {}
                for nb in names_all:
                    pb = compb[nb]
                    if pb:
                        Ecache[nb] = fr.einv(fr.t(a, pb))
                        smA[nb] = fr.Pi(fr.sm(sa_, la, pb))
                for x, na in enumerate(names_all):
                    pa = compa[na]
                    if not pa:
                        continue
                    Tb_pa = fr.t(bb, pa)
                    smb_pa = fr.sm(sb_, lb, pa)
                    for y, nb in enumerate(names_all):
                        if nb not in Ecache:
                            continue
                        Q0[x][y] += sz * coef * sab * fr.dot(smb_pa, smA[nb])
                        addE(x, y, [(sk, -sz * b * sab * fr.dot(Tb_pa, vv)) for sk, vv in Ecache[nb]])
    # vacuum-vacuum
    if C.mode == 'sym':
        fr = Frame(C, ['M'])
        ph = {nm: fr.phi(nm) for nm in names_all}
        for x, na in enumerate(names_all):
            for y, nb in enumerate(names_all):
                if ph[na] and ph[nb]:
                    if not C.sym:
                        G[x][y] += fr.dot(ph[na], ph[nb])
                    Q0[x][y] += 2 * fr.dot(fr.Pi(ph[na]), ph[nb])
                    Ua, Ub = fr.jp(ph[na]), fr.jp(ph[nb])
                    if Ua and Ub:
                        addE(x, y, [(sk, -b * fr.dot(Ua, vv)) for sk, vv in fr.einv(Ub)])
    keep, norms = [], []
    for x, nm in enumerate(names_all):
        nrm = norm_of(C, nm)
        if nrm != 0:
            keep.append(x)
            norms.append(nrm)
    if C.sym:
        Q = []
        for x in keep:
            row = []
            for y in keep:
                val = RF(Q0[x][y]) if not isinstance(Q0[x][y], int) else RF(fq(Q0[x][y]))
                for sk, M in QE.items():
                    if M[x][y] != 0:
                        val = val + RF(M[x][y], C.D[sk])
                row.append(val)
            Q.append(row)
    else:
        Q = []
        for x in keep:
            row = []
            for y in keep:
                val = Q0[x][y]
                for sk, M in QE.items():
                    val = val + M[x][y]
                row.append(val)
            Q.append(row)
    return [names_all[x] for x in keep], norms, Q


def comb_ok(n, r):
    if n < 0 or r < 0 or r > n:
        return 0
    return comb(n, r)


def norm_of(C, nm):
    """squared norm of the twisted vector nm (number of states in its support, all coefficients +-1)"""
    i, m, q0, tp = C.i, C.m, C.q0, C.taup
    p2 = 2 ** i
    if C.case == 'MT':
        if C.mode == 'sym':
            if nm in ('A0', 'A1'):
                x = int(nm[1])
                return m * p2 * comb_ok(m - 1, q0 - x) if m >= 1 else 0
            if nm[0] == 'B':
                x, y = int(nm[1]), int(nm[2])
                return 2 * i * (p2 // 2) * comb_ok(m, q0 + 1 - x - y) if i >= 1 else 0
            if nm == 'C':
                return tp * p2 * comb_ok(m, q0)
            if nm in ('D0', 'D1'):
                z = int(nm[1])
                return tp * p2 * comb_ok(m, q0 + 1 - z)
            if nm == 'CM':
                return p2 * comb_ok(m, q0)
            if nm in ('DM0', 'DM1'):
                z = int(nm[2])
                return p2 * comb_ok(m, q0 + 1 - z)
            if nm == 'G':
                return p2 * comb_ok(m, q0 + 1) if C.j + 1 <= C.k else 0
        else:
            if nm == 'Ca':
                return 2 * p2 * comb_ok(m, q0)
            z = int(nm[2])
            return 2 * p2 * comb_ok(m, q0 + 1 - z)
    else:
        if C.mode == 'sym':
            if nm[0] == 'O':
                y0 = int(nm[1])
                return p2 * comb_ok(m, q0 - y0)
            if nm[0] == 'A':
                x, y0 = int(nm[1]), int(nm[2])
                return m * p2 * comb_ok(m - 1, q0 - x - y0) if m >= 1 else 0
            if nm[0] == 'B':
                x, y, y0 = int(nm[1]), int(nm[2]), int(nm[3])
                return 2 * i * (p2 // 2) * comb_ok(m, q0 + 1 - x - y - y0) if i >= 1 else 0
            if nm[0] == 'C':
                y0 = int(nm[1])
                return tp * p2 * comb_ok(m, q0 - y0)
            if nm[0] == 'D':
                z, y0 = int(nm[1]), int(nm[2])
                return tp * p2 * comb_ok(m, q0 + 1 - z - y0)
            if nm[0] == 'G':
                y0 = int(nm[1])
                return p2 * comb_ok(m, q0 + 1 - y0) if C.j + 1 <= C.k else 0
        else:
            if nm[:2] == 'Ca':
                y0 = int(nm[2])
                return 2 * p2 * comb_ok(m, q0 - y0)
            z, y0 = int(nm[2]), int(nm[3])
            return 2 * p2 * comb_ok(m, q0 + 1 - z - y0)
    raise ValueError(nm)


def kernel_vectors(C, names):
    """twisted parts of the known kernel p (x) H_j of the block (rescaled frame): list of (label, {name: coef})"""
    vecs = []
    if C.mode == 'anti':
        return vecs
    if C.case == 'MT':
        v = {}
        for nm, cf in (('A0', 1), ('A1', 1), ('B10', 1), ('B01', -1), ('C', 1), ('CM', 1)):
            if nm in names:
                v[nm] = cf
        if v:
            vecs.append(('p x h', v))
        return vecs
    for y0 in (0, 1):
        v = {}
        for nm, cf in (('O%d' % y0, 1), ('A0%d' % y0, 1), ('A1%d' % y0, 1), ('B10%d' % y0, 1), ('B01%d' % y0, -1),
                       ('C%d' % y0, 1)):
            if nm in names:
                v[nm] = cf
        if v:
            vecs.append(('p x h%d' % y0, v))
    return vecs


def i_range(k, tau, case):
    up = k - tau if case == 'MT' else k - tau - 1
    return list(range(0, max(up, 0) // 2 + 1))


def all_forms(k, jmin=3, jmax=None):
    """every (j, tau, case, mode, i) of a type I block with a letter 1, 3 <= j <= k - 1 (or the given range)"""
    if jmax is None:
        jmax = k - 1
    out = []
    for j in range(jmin, jmax + 1):
        for tau in range(1, j + 1):
            for case in ('MT', 'MU'):
                if case == 'MU' and k - tau < 1:
                    continue
                if case == 'MT' and tau > k:
                    continue
                for mode in ('sym', 'anti'):
                    for i in i_range(k, tau, case):
                        C = Calc(k, j, tau, i, case, mode, sigma=1)
                        if not C.valid():
                            continue
                        if all(norm_of(C, nm) == 0 for nm in NAMES[(case, mode)]):
                            continue
                        out.append((j, tau, case, mode, i))
    return out
