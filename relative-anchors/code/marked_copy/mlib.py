"""Orbit calculus for the letter-1 block of (R_inf)_j with one marked copy.

Copy 0 (label 'M') has squared weight sigma > 0, the u = k - 1 copies of U = {1, ..., u} have squared weight 1.
Everything is in the rescaled frame of ylib.py (each state divided by prod_t w_t^{N_t}), in units 2^Z of Lemma B,
where (README of relative-anchors, Lemma A of code/middle_j2_l1 adapted to one marked copy):
  Y~[(a,y'),(b,y)] = (1/2 - [a = b]) s_a s_b pi_z [y' = y - a + b] (z = y - a) + [a = b][y' = y] s_a pi_y (2 - n_a)
                     - s_a s_b <y'| t_b^T E^^{-1} t_a |y>,
  Y~[(a,y'),(vac,y)] = -[y' = y - a] pi_y + s_a <y'| jb^- E^^{-1} t_a |y>,
  Y~[(vac,y'),(vac,y)] = 2 [y' = y] pi_y - <y'| jb^- E^^{-1} jb^+ |y>,
with s_M = sigma, s_t = 1 on U, pi_y = sigma^{y_0}, jb^{+-} the unweighted raising/lowering operators,
t_a = 1/2 sigma_a^- jb^+ + N_a (N_a = 3/2 |0><0| + 1/2 |1><1|) and E^ = D_w^{-1} E D_w^{-1}.  On a level c split
into A (letter 0 at M; c-subsets of U) and B (letter 1 at M; (c-1)-subsets of U):
  E^ = [[2 sigma + 2(u - c) + 1/2 J_U^+ J_U^-,   1/2 J_U^+ sigma_M^-],
        [1/2 sigma_M^+ J_U^-,   (2(u - c + 1) + sigma/2 + 1/2 J_U^+ J_U^-)/sigma]].
On the S_u-irrep (u - i', i') this is a 2 x 2 matrix; its inverse is applied by Newton interpolation in DU = J_U^- J_U^+
over the nodes i' of a window around i (the run asserts that the window is exact: the remainder vanishes).

Symmetry: H_i = (S_2 wr S_i) x S_m on U (pairs (1,2), (3,4), ...; m = u - 2i), character chi_i = sign on each pair
swap.  Twisted invariants (rescaled coordinates, coefficients +-1), with q0 = j - i:
  O{y0}:        row 1@M,  column y0 at M, singlets on all pairs, Dicke(R, q0 - y0);
  A{x}{y0}:     row 1@a, a in R, column x at a, y0 at M, singlets, Dicke(R - a, q0 - y0 - x);
  B{x}{y}{y0}:  row 1@a, a in a pair (a, a'), column x at a, y at a', y0 at M, singlets on the other pairs,
                Dicke(R, q0 + 1 - x - y - y0), sign -1 when a is the second site of its pair;
  G{z0}:        vacuum row, column z0 at M, singlets, Dicke(R, q0 + 1 - z0).
A vector in a frame (distinguished sites: M, up to two U-sites with their partners) is a dict
(letters on the frame, d) -> coefficient, d the offset of the Dicke count on R - frame from q0.

Modes: Calc(k, j, i, sigma=Fraction) is exact numeric (python-flint fmpq); Calc(k, j, i, sigma=None) is symbolic in
sigma: coefficients in Q[sigma] (fmpq_poly), every E^-inverse at level c carries the denominator D_c (a product of the
2 x 2 determinants), and the forms are returned as reduced rational functions (ufield.RF).
bscale multiplies every term that contains an inverse E^-1 (the B-part), as in j1lib.py.
"""
import itertools
from fractions import Fraction as Fr
from math import comb
import flint
from ufield import RF

fq = flint.fmpq
PL = flint.fmpq_poly

NAMES = ['O0', 'O1', 'A00', 'A01', 'A10', 'A11',
         'B000', 'B001', 'B010', 'B011', 'B100', 'B101', 'B110', 'B111', 'G0', 'G1']


def tofq(x):
    if isinstance(x, Fr):
        return fq(x.numerator, x.denominator)
    return fq(x)


class Calc:
    def __init__(self, k, j, i, sigma=None, bscale=1, full_nodes=False):
        self.k, self.j, self.i = k, j, i
        self.u = k - 1
        self.m = self.u - 2 * i
        self.q0 = j - i
        self.sym = sigma is None
        if self.sym:
            self.sig = PL([0, 1])
            self.one = PL([1])
        else:
            self.sig = tofq(sigma)
            self.one = fq(1)
        self.zero = self.one * 0
        self.bscale = tofq(bscale)
        self.full_nodes = full_nodes
        self._einv = {}
        self.window_failures = 0

    # ---------------- E^-inverse data per level ----------------
    def lam(self, ell, ip):
        """eigenvalue of DU = J_U^- J_U^+ on the irrep (u - i', i') of the ell-subsets of U"""
        return (self.u - ell - ip) * (ell + 1 - ip)

    def irr(self, ell):
        """irreps i' present on the ell-subsets of U"""
        if ell < 0 or ell > self.u:
            return []
        return list(range(0, min(ell, self.u - ell) + 1))

    def order(self, nodes):
        w = [x for x in nodes if self.i - 2 <= x <= self.i + 2]
        if self.full_nodes:
            return w + [x for x in nodes if x not in w]
        return w

    def einv_data(self, c):
        if c in self._einv:
            return self._einv[c]
        u, s = self.u, self.sig
        IA = self.irr(c)          # A: U-level c
        IB = self.irr(c - 1)      # B: U-level c - 1
        vals = {}
        for ip in sorted(set(IA) | set(IB)):
            mu = (u - c + 1 - ip) * (c - ip)
            nu = (u - c + 2 - ip) * (c - 1 - ip)
            a = 2 * s + tofq(Fr(2 * (u - c)) + Fr(mu, 2))
            e = s * fq(1, 2) + tofq(Fr(2 * (u - c + 1)) + Fr(nu, 2))
            if ip in IA and ip in IB:
                D = a * e - s * fq(mu, 4)
                vals[ip] = ('both', a, e, D, mu)
            elif ip in IA:
                vals[ip] = ('A', a, e, a, mu)
            else:
                vals[ip] = ('B', a, e, e, mu)
        nA = self.order(IA)
        nB = self.order(IB)
        used = sorted(set(nA) | set(nB))
        if self.sym:
            Dc = self.one
            for ip in used:
                Dc = Dc * vals[ip][3]
        else:
            Dc = self.one

        def fval(kind, ip):
            t, a, e, D, mu = vals[ip]
            if kind == 1:
                num, den = (e, D) if t == 'both' else (self.one, a)
            elif kind == 2:
                num, den = ((-s * fq(1, 2)), D) if t == 'both' else (self.zero, self.one)
            else:
                num, den = (s * a, D) if t == 'both' else (s, e)
            return num, den

        def newton(kind, nodes, ell):
            lamv = [self.lam(ell, ip) for ip in nodes]
            coefs = []
            diag = []
            for n, ip in enumerate(nodes):
                num, den = fval(kind, ip)
                if self.sym:
                    fv = RF(num, den)
                else:
                    fv = num / den
                new = [fv]
                for l in range(1, n + 1):
                    new.append((new[l - 1] - diag[l - 1]) / (lamv[n] - lamv[n - l]))
                diag = new
                coefs.append(new[n])
            if self.sym:
                out = []
                for cf in coefs:
                    x = cf * RF(Dc)
                    assert x.d == PL([1]), 'Newton coefficient times D_c is not a polynomial'
                    out.append(x.n)
                coefs = out
            return lamv, coefs

        data = {'A1': newton(1, nA, c), 'B2': newton(2, nB, c - 1), 'B3': newton(3, nB, c - 1), 'D': Dc}
        self._einv[c] = data
        return data


class Frame:
    def __init__(self, calc, labels):
        """labels: 'M' first, then U-site labels like 'P0a', 'P0b', 'R0' (pair labels come with partners)"""
        self.C = calc
        self.labels = labels
        assert labels[0] == 'M'
        self.pairs = sorted({l[1] for l in labels if l[0] == 'P'})
        self.rs = [l for l in labels if l[0] == 'R']
        self.pD = len(self.pairs)
        self.nprime = calc.m - len(self.rs)
        self.undist = calc.i - self.pD
        self.pos = {l: n for n, l in enumerate(labels)}
        self.upos = [n for n, l in enumerate(labels) if l != 'M']

    # ---- helpers
    def inrange(self, d):
        q = self.C.q0 + d
        return 0 <= q <= self.nprime

    def level(self, key):
        L, d = key
        return sum(L) + self.C.q0 + d + self.undist

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

    # ---- operators (unweighted)
    def jUp(self, v):
        out = {}
        q0 = self.C.q0
        for (L, d), c in v.items():
            for s in self.upos:
                if L[s] == 0:
                    self.addto(out, (L[:s] + (1,) + L[s + 1:], d), c)
            if self.inrange(d + 1):
                self.addto(out, (L, d + 1), c * (q0 + d + 1))
        return out

    def jUm(self, v):
        out = {}
        q0 = self.C.q0
        for (L, d), c in v.items():
            for s in self.upos:
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

    def sm(self, lab, v):
        s = self.pos[lab]
        return {(L[:s] + (0,) + L[s + 1:], d): c for (L, d), c in v.items() if L[s] == 1}

    def nop(self, lab, v):
        s = self.pos[lab]
        return {key: c for key, c in v.items() if key[0][s] == 1}

    def Nop(self, lab, v):
        s = self.pos[lab]
        return {(L, d): c * (fq(1, 2) if L[s] == 1 else fq(3, 2)) for (L, d), c in v.items()}

    def t(self, lab, v):
        out = self.scale(self.sm(lab, self.jp(v)), fq(1, 2))
        return self.add(out, self.Nop(lab, v))

    def Pi(self, v):
        s = self.C.sig
        return {(L, d): (c * s if L[0] == 1 else c) for (L, d), c in v.items()}

    def DU(self, v):
        return self.jUm(self.jUp(v))

    def einv(self, v):
        """E^-inverse of v (all keys at one level c); in symbolic mode the result is D_c times the inverse."""
        if not v:
            return {}
        lv = {self.level(key) for key in v}
        assert len(lv) == 1, lv
        c = lv.pop()
        data = self.C.einv_data(c)
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
            self.C.window_failures += 1
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
                d = key[1]
                wgt = (2 ** und) * comb(self.nprime, q0 + d)
                s = s + x * y * wgt
        return s

    # ---- basis vectors
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
                key = (tuple(L), d)
                self.addto(out, key, self.C.one * coef)
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
    """H_i-orbits of sites and of ordered pairs of sites of {M} u U, with their sizes"""
    m, i = C.m, C.i
    sr = [('M', 1)]
    if m >= 1:
        sr.append(('R0', m))
    if i >= 1:
        sr.append(('P0a', 2 * i))
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
    return sr, pr


def reduced(C):
    """the reduced form on the nonzero twisted vectors: returns (names, norms, Q).
    numeric mode: Q entries fmpq; symbolic mode: Q entries RF in sigma; norms are the squared norms (fmpq)."""
    n = len(NAMES)
    sr, pr = orbits(C)
    j, k = C.j, C.k
    sig = C.sig
    # accumulators: base polynomial part and parts over D_j, D_{j+1}, D_{j+2}
    Q0 = [[C.zero] * n for _ in range(n)]
    QE = {c: [[C.zero] * n for _ in range(n)] for c in (j, j + 1, j + 2)}
    G = [[fq(0)] * n for _ in range(n)]
    b = C.bscale
    s_of = lambda lab: sig if lab == 'M' else C.one
    # site terms
    for (a, sz) in sr:
        fr = Frame(C, frame_labels(a))
        comps = {nm: fr.psi(nm, a) for nm in NAMES}
        phis = {nm: fr.phi(nm) for nm in NAMES}
        sa = s_of(a)
        Cbase, CE = {}, {}
        for nm in NAMES:
            ph = phis[nm]
            if ph:
                Cbase[nm] = fr.scale(fr.sm(a, fr.Pi(ph)), -1)
                CE[nm] = fr.scale(fr.jm(fr.einv(fr.t(a, ph))), sa * b)
        for x, na in enumerate(NAMES):
            pa = comps[na]
            if not pa:
                continue
            Pipa = fr.Pi(pa)
            for y, nb in enumerate(NAMES):
                pb = comps[nb]
                if pb:
                    G[x][y] += sz * (fr.dot(pa, pb) if not C.sym else 0)
                    Q0[x][y] += sz * sa * (2 * fr.dot(Pipa, pb) - fr.dot(Pipa, fr.nop(a, pb)))
                if nb in Cbase:
                    Q0[x][y] += sz * fr.dot(pa, Cbase[nb])
                    QE[j + 1][x][y] += sz * fr.dot(pa, CE[nb])
                    Q0[y][x] += sz * fr.dot(pa, Cbase[nb])
                    QE[j + 1][y][x] += sz * fr.dot(pa, CE[nb])
    # pair terms
    for (a, bb, sz) in pr:
        fr = Frame(C, frame_labels(a, bb))
        compa = {nm: fr.psi(nm, a) for nm in NAMES}
        compb = {nm: fr.psi(nm, bb) for nm in NAMES}
        sab = s_of(a) * s_of(bb)
        coef = fq(1, 2) - (1 if a == bb else 0)
        Ecache = {}
        smA = {}
        for nb in NAMES:
            pb = compb[nb]
            if pb:
                Ecache[nb] = fr.einv(fr.t(a, pb))
                smA[nb] = fr.Pi(fr.sm(a, pb))
        for x, na in enumerate(NAMES):
            pa = compa[na]
            if not pa:
                continue
            Tb_pa = fr.t(bb, pa)
            smb_pa = fr.sm(bb, pa)
            for y, nb in enumerate(NAMES):
                if nb not in Ecache:
                    continue
                Q0[x][y] += sz * coef * sab * fr.dot(smb_pa, smA[nb])
                QE[j][x][y] += -sz * b * sab * fr.dot(Tb_pa, Ecache[nb])
    # vacuum-vacuum
    fr = Frame(C, ['M'])
    ph = {nm: fr.phi(nm) for nm in NAMES}
    for x, na in enumerate(NAMES):
        for y, nb in enumerate(NAMES):
            if ph[na] and ph[nb]:
                if not C.sym:
                    G[x][y] += fr.dot(ph[na], ph[nb])
                Q0[x][y] += 2 * fr.dot(fr.Pi(ph[na]), ph[nb])
                Ua, Ub = fr.jp(ph[na]), fr.jp(ph[nb])
                if Ua and Ub:
                    QE[j + 2][x][y] += -b * fr.dot(Ua, fr.einv(Ub))
    # norms (independent of sigma): recompute in numeric terms
    keep = []
    norms = []
    for x, nm in enumerate(NAMES):
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
                for c in (j, j + 1, j + 2):
                    if QE[c][x][y] != 0:
                        val = val + RF(QE[c][x][y], C.einv_data(c)['D'])
                row.append(val)
            Q.append(row)
    else:
        Q = [[Q0[x][y] + QE[j][x][y] + QE[j + 1][x][y] + QE[j + 2][x][y] for y in keep] for x in keep]
    return [NAMES[x] for x in keep], norms, Q


def norm_of(C, nm):
    """squared norm of the twisted vector nm (number of states times 1), by direct counting"""
    i, m, q0 = C.i, C.m, C.q0
    if nm[0] == 'O':
        y0 = int(nm[1])
        return (2 ** i) * comb_ok(m, q0 - y0)
    if nm[0] == 'A':
        x, y0 = int(nm[1]), int(nm[2])
        return m * (2 ** i) * comb_ok(m - 1, q0 - y0 - x) if m >= 1 else 0
    if nm[0] == 'B':
        x, y, y0 = int(nm[1]), int(nm[2]), int(nm[3])
        return 2 * i * (2 ** (i - 1)) * comb_ok(m, q0 + 1 - x - y - y0) if i >= 1 else 0
    z0 = int(nm[1])
    if C.j + 1 > C.k:
        return 0
    return (2 ** i) * comb_ok(m, q0 + 1 - z0)


def comb_ok(n, r):
    if n < 0 or r < 0 or r > n:
        return 0
    return comb(n, r)


# ---------------- known kernel vectors ----------------
def kernel_vectors(names, k, j, i, sig):
    """the twisted invariants of the known kernel: list of (label, dict name -> coefficient); sig a field element"""
    u = k - 1
    vecs = []
    for y0 in (0, 1):
        v = {}
        for nm, cf in (('O%d' % y0, 1), ('A0%d' % y0, 1), ('A1%d' % y0, 1), ('B10%d' % y0, 1), ('B01%d' % y0, -1)):
            if nm in names:
                v[nm] = cf
        if v:
            vecs.append(('p x h%d' % y0, v))
    if j == k - 2 and i == 0:
        v = {nm: cf for nm, cf in (('G0', sig), ('G1', 1)) if nm in names}
        vecs.append(('K2', v))
    if j == k - 1 and i == 0:
        v = {}
        for nm, cf in (('O0', u * sig), ('O1', u), ('A10', -sig * sig), ('A01', -sig), ('A11', -sig)):
            if nm in names:
                v[nm] = cf
        vecs.append(('T0', v))
    if j == k - 1 and i == 1:
        v = {}
        for nm, cf in (('B110', sig), ('B011', 1), ('B101', 1), ('B111', 1)):
            if nm in names:
                v[nm] = cf
        vecs.append(('T1', v))
    return vecs
