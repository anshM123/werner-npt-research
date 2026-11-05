"""Family 2: the type I blocks of (R_inf)_2 with T = {t} (two letters 1), the U-weights on two levels, exact and
symbolic in the level sizes (pattern sums as in middle_j1_l1/l1sym.py).

Copies: t (squared weight a; level 2 below, of size 1) and U = U_0 + U_1, |U_i| = u_i, squared weight Q_i on U_i.
L = a + u_0 Q_0 + u_1 Q_1.  Units 2^{k-2}; rescaled frame (each state multiplied by the product of the weights over its
letters 1).  States:
  R_P     = (2@t, 1@P),      P a pair of [k];
  P_{r,c} = (1@r, 2@t 1@c),  r in [k], c in U;
  Q_P     = (vac, 2@t 1@P),  P a pair of U.
Closed forms (NOTES.md, Section 3; Delta_x = 2L - a - 2 s_x, gamma = sum_{x in U} s_x/Delta_x, sg = 2 + gamma,
nu_{r,c} = 2 if r not in {t, c}, 1 if r = t, 0 if r = c; C(P) = 2(L - s_P), C3(T) = 2(L - s_T), C'(P) = 2L - a - 2 s_P,
C4(T) = 2L - a - 2 s_T; eps_x = -1 if x = t, else 1):
  RR: diag s_P (2 - [t in P]);
  PP: (1/2 - [r=r'])[c'=r][c=r'] s_r s_r' + [r=r'][c=c'](2 - [r in {t,c}]) s_r s_c - b beta,
      beta = [r'=c][r=c'] s_c s_c' gamma/(2 sg) + [r'=c] s_r s_c s_c' nu_{r,c'}/(Delta_c' sg)
             + [r=c'] s_c' s_r' s_c nu_{r',c}/(Delta_c sg) + s_r s_r' nu_{r',c} nu_{r,c'} ([c=c'] s_c/Delta_c - s_c s_c'/(Delta_c Delta_c' sg));
  RP: (1/2 - [r=t])[P = {c,r}] s_r s_c - b ( 1/4 [r in P] s_r Mt[P-r, c]
      + [r not in P] s_r/C(P) ([c in P] eps_{P-c} s_{P-c} s_c - 1/2 (s_z Mt[x,c] + s_x Mt[z,c])) ),  P = {x, z},
      Mt = (A^-1 N_{0,t-})~ on the singles of [k], A = I + N/2, N = J^- C^-1 J^+ (pairs C), N_{0,t-} with w_t -> -w_t
      in J^+;
  PQ: -[r in P] s_r At[c, P-r] + [r not in P] s_r nu_{r,P} (s_z At[c,x] + s_x At[c,z])/C'(P)   (b = 1),
      At = (A_t^-1)~ on U, A_t = I + N'/2, N' = J^- C'^-1 J^+ on the singles of the t-sector; A-part
      -[r in P][c = P - r] s_r s_c;
  RQ: -[P = P'] s_P + 1/2 (A2^-1 N2_{0,t-})~[P, P'],  A2 = I + N2/2 on the pairs of [k], N2 = J^- C3^-1 J^+;
  QQ: 2 W4^-1, W4 = D^-1 (I + N4/2) D^-1 on the pairs of U, N4 = J^- C4^-1 J^+ (t-sector).
Labels: (level, name, id); level 2 is the copy t (always named); anonymous labels only on levels 0 and 1.
"""
from fractions import Fraction as Fr
from l1sym import num, ff, canon, from_key, mm, mT, minv, det


class Model2:
    def __init__(self, one, a, Q, U, b=1):
        self.one, self.zero = one, one * 0
        self.Q = [Q[0], Q[1], a]
        self.U = [U[0], U[1], one]
        self.a = a
        self.b = num(one, b)
        z = self.zero
        L = U[0] * Q[0] + U[1] * Q[1] + a
        self.L = L
        sv = self.Q
        Uv = self.U
        # t-sector singles (U): Delta, gamma, sigma
        self.gam = sum((Uv[i] * sv[i] / (2 * L - a - 2 * sv[i]) for i in (0, 1)), z)
        self.sg = self.gam + 2
        # ---- small invariant algebra on the singles of [k] (orbits 2 = t, 0, 1) and of U (orbits 0, 1)
        self.Ach = self._inv_singles(levels=(2, 0, 1), C=lambda x, y: 2 * (L - x - y), diag_extra=None)
        self.Atch = self._inv_singles(levels=(0, 1), C=lambda x, y: 2 * L - a - 2 * x - 2 * y, diag_extra=None)
        # Mt = Ach * Ncheck_{0,t-} * D^2  on the singles of [k]
        self.Mt = self._Mt()

    # ---------------- helpers on label-valued invariant operators of a single-index space
    def _size(self, lv):
        return self.U[lv]

    def _sumexcl(self, levels, excl, fun):
        tot = self.zero
        for lv in levels:
            n = self.U[lv] - (1 if lv == excl else 0)
            if n != 0:
                tot = tot + n * fun(lv)
        return tot

    def _to_red(self, levels, f):
        """reduced matrix and std scalars of the invariant op with label values f(ly, lx, same)."""
        M = []
        for ly in levels:
            row = []
            for lx in levels:
                if ly == lx and lx == 2:
                    row.append(f(2, 2, True))
                elif ly == lx:
                    # a level of fixed size 0 or 1 has no pair of distinct copies: that value is not used
                    off = self.zero if (self.U[ly] - 1 == 0 or self.U[ly] == 0) else (self.U[ly] - 1) * f(ly, ly, False)
                    row.append(f(ly, ly, True) + off)
                else:
                    row.append(self.U[lx] * f(ly, lx, False))
            M.append(row)
        sd = {lv: f(lv, lv, True) - (f(lv, lv, False) if (self.U[lv] - 1 != 0 and self.U[lv] != 0) else self.zero)
              for lv in levels if lv != 2}
        return M, sd

    def _from_red(self, levels, M, sd):
        """label values of the invariant op with reduced matrix M and std scalars sd."""
        vals = {}
        for i, ly in enumerate(levels):
            for j, lx in enumerate(levels):
                if ly == lx and lx == 2:
                    vals[(2, 2, True)] = M[i][j]
                elif ly == lx:
                    u = self.U[ly]
                    if u == 0:
                        vals[(ly, ly, True)] = vals[(ly, ly, False)] = self.zero
                        continue
                    ep = (M[i][j] - sd[ly]) / u
                    vals[(ly, ly, True)] = sd[ly] + ep
                    vals[(ly, ly, False)] = ep
                else:
                    vals[(ly, lx, False)] = self.zero if self.U[lx] == 0 else M[i][j] / self.U[lx]
        return vals

    def _inv_singles(self, levels, C, diag_extra):
        """rescaled inverse (D A^-1 D) of A = I + N/2 on the singles of the given levels, N = J^- C^-1 J^+ with
        C(s_x, s_y) the energy of the pair {x, y}: Ahat = D^-1 A D^-1 has diagonal (1 + N[x,x]/2)/s_x and
        off-diagonal 1/(2 C)."""
        s = self.Q

        def f(ly, lx, same):
            if same:
                nxx = self._sumexcl(levels, ly, lambda lv: s[lv] / C(s[ly], s[lv]))
                return (self.one + nxx / 2) / s[ly]
            return self.one / (2 * C(s[ly], s[lx]))
        M, sd = self._to_red(levels, f)
        Mi = minv(M, self.one)
        sdi = {lv: self.one / v for lv, v in sd.items()}
        return self._from_red(levels, Mi, sdi)

    def _Mt(self):
        levels = (2, 0, 1)
        s = self.Q
        L = self.L
        C = lambda x, y: 2 * (L - x - y)
        eps = lambda lv: -1 if lv == 2 else 1

        def ncheck(ly, lx, same):          # Ncheck_{0,t-}[y, c] = N_{0,t-}[y, c]/(w_y w_c)
            if same:
                return self._sumexcl(levels, ly, lambda lv: eps(lv) * s[lv] / C(s[ly], s[lv])) / s[ly]
            return eps(ly) * self.one / C(s[ly], s[lx])
        Ma, sda = self._to_red(levels, lambda ly, lx, same: self.Ach[(ly, lx, same)] if (ly, lx, same) in self.Ach
                               else self.zero)
        Mn, sdn = self._to_red(levels, ncheck)
        Md, sdd = self._to_red(levels, lambda ly, lx, same: s[ly] if same else self.zero)
        M = mm(mm(Ma, Mn, self.zero), Md, self.zero)
        sd = {lv: sda[lv] * sdn[lv] * sdd[lv] for lv in sda}
        return self._from_red(levels, M, sd)

    # ---------------- entries on labels
    def s(self, x):
        return self.Q[x[0]]

    def lv(self, x, y):
        """the (ly, lx, same) key of the pair of labels."""
        return (x[0], y[0], x == y)

    def RR(self, P1, P2):
        x, z = P1
        if {x, z} != set(P2):
            return self.zero
        return self.s(x) * self.s(z) * (2 - (1 if (x[0] == 2 or z[0] == 2) else 0))

    def nu(self, r, c):
        if r == c:
            return 0
        return 1 if r[0] == 2 else 2

    def PP(self, r, c, rp, cp):
        s = self.s
        A = self.zero
        if cp == r and c == rp:
            A = A + (self.one / 2 - (1 if r == rp else 0)) * s(r) * s(rp)
        if r == rp and c == cp:
            A = A + (2 - (1 if (r[0] == 2 or r == c) else 0)) * s(r) * s(c)
        D = lambda x: 2 * self.L - self.a - 2 * s(x)
        g, sg = self.gam, self.sg
        be = self.zero
        if rp == c and r == cp:
            be = be + s(c) * s(cp) * g / (2 * sg)
        if rp == c:
            be = be + s(r) * s(c) * s(cp) * self.nu(r, cp) / (D(cp) * sg)
        if r == cp:
            be = be + s(cp) * s(rp) * s(c) * self.nu(rp, c) / (D(c) * sg)
        n1, n2 = self.nu(rp, c), self.nu(r, cp)
        if n1 and n2:
            inner = -s(c) * s(cp) / (D(c) * D(cp) * sg)
            if c == cp:
                inner = inner + s(c) / D(c)
            be = be + s(r) * s(rp) * n1 * n2 * inner
        return A - self.b * be

    def Mtv(self, x, c):
        return self.Mt[self.lv(x, c)]

    def RP(self, x, z, r, c):
        s = self.s
        A = self.zero
        if r != c and {x, z} == {c, r}:
            A = (self.one / 2 - (1 if r[0] == 2 else 0)) * s(r) * s(c)
        B = self.zero
        if r == x or r == z:
            o = z if r == x else x
            B = B + s(r) * self.Mtv(o, c) / 4
        else:
            CP = 2 * (self.L - s(x) - s(z))
            inner = -(s(z) * self.Mtv(x, c) + s(x) * self.Mtv(z, c)) / 2
            if c == x or c == z:
                o = z if c == x else x
                inner = inner + (-1 if o[0] == 2 else 1) * s(o) * s(c)
            B = B + s(r) * inner / CP
        return A - self.b * B

    def Atv(self, c, x):
        return self.Atch[self.lv(c, x)]

    def PQ(self, r, c, x, z):
        s = self.s
        if r == x or r == z:
            o = z if r == x else x
            full = -s(r) * self.Atv(c, o)
            Ap = -s(r) * s(c) if c == o else self.zero
        else:
            nu = 1 if r[0] == 2 else 2
            full = s(r) * nu * (s(z) * self.Atv(c, x) + s(x) * self.Atv(c, z)) / (2 * self.L - self.a - 2 * s(x) - 2 * s(z))
            Ap = self.zero
        return Ap + self.b * (full - Ap)

    # pair-space operators (for RQ and QQ, through reduced forms)
    def A2hat(self, x, z, x2, z2):
        """Ahat2 = D^-1 (I + N2/2) D^-1 on the pairs of [k]."""
        s = self.s
        P, P2 = {x, z}, {x2, z2}
        L = self.L
        if P == P2:
            n2 = self.zero
            for lv in (0, 1, 2):
                cnt = sum(1 for y in P if y[0] == lv)
                if self.U[lv] - cnt != 0:
                    n2 = n2 + (self.U[lv] - cnt) * self.Q[lv] / (2 * (L - s(x) - s(z) - self.Q[lv]))
            return (self.one + n2 / 2) / (s(x) * s(z))
        I = P & P2
        if len(I) == 1:
            c = next(iter(I))
            y1 = next(iter(P - I))
            y2 = next(iter(P2 - I))
            return self.one / (4 * s(c) * (L - s(y1) - s(y2) - s(c)))
        return self.zero

    def N2c(self, x, z, x2, z2):
        """Ncheck2_{0,t-}[P, P'] = N2_{0,t-}[P, P']/(w_P w_P'), P = {x, z} a pair of [k], P' = {x2, z2} a pair of U."""
        s = self.s
        P, P2 = {x, z}, {x2, z2}
        L = self.L
        if P == P2:
            tot = self.zero
            for lv in (0, 1, 2):
                cnt = sum(1 for y in P if y[0] == lv)
                if self.U[lv] - cnt != 0:
                    tot = tot + (self.U[lv] - cnt) * (-1 if lv == 2 else 1) * self.Q[lv] / (2 * (L - s(x) - s(z) - self.Q[lv]))
            return tot / (s(x) * s(z))
        I = P & P2
        if len(I) == 1:
            c = next(iter(I))
            y1 = next(iter(P - I))       # the element of P not in P'
            y2 = next(iter(P2 - I))
            return (-1 if y1[0] == 2 else 1) * self.one / (2 * s(c) * (L - s(y1) - s(y2) - s(c)))
        return self.zero

    def W4(self, x, z, x2, z2):
        """W4 = D^-1 (I + N4/2) D^-1 on the pairs of U (t-sector, C4(T) = 2L - a - 2 s_T)."""
        s = self.s
        P, P2 = {x, z}, {x2, z2}
        L, a = self.L, self.a
        if P == P2:
            n2 = self.zero
            for lv in (0, 1):
                cnt = sum(1 for y in P if y[0] == lv)
                if self.U[lv] - cnt != 0:
                    n2 = n2 + (self.U[lv] - cnt) * self.Q[lv] / (2 * L - a - 2 * (s(x) + s(z) + self.Q[lv]))
            return (self.one + n2 / 2) / (s(x) * s(z))
        I = P & P2
        if len(I) == 1:
            c = next(iter(I))
            y1 = next(iter(P - I))
            y2 = next(iter(P2 - I))
            return self.one / (2 * s(c) * (2 * L - a - 2 * (s(y1) + s(y2) + s(c))))
        return self.zero


# ---------------------------------------------------------------- test vectors
T_LAB = (2, 't', 't')


def nm(name, lev):
    return (lev, name, name)


def fo(x, plus, minus):
    return 1 if x[1] == plus else (-1 if x[1] == minus else 0)


def part_spec2(part):
    """(named labels, R-vectors, P-vectors, Q-vectors, kernel coefficient lists on the P-vectors, gauge P-indices,
    names).  The irrep sits on level 0 ('mix': both).  R- and Q-vectors are functions of an ordered pair (x, z)
    (symmetric); P-vectors of (r, c)."""
    isU = lambda x: x[0] in (0, 1)
    if part == 'inv':
        named = [T_LAB]
        lv = lambda x, i: x[0] == i
        R = [lambda x, z: int((lv(x, 2) and lv(z, 0)) or (lv(z, 2) and lv(x, 0))),
             lambda x, z: int((lv(x, 2) and lv(z, 1)) or (lv(z, 2) and lv(x, 1))),
             lambda x, z: int(lv(x, 0) and lv(z, 0)),
             lambda x, z: int((lv(x, 0) and lv(z, 1)) or (lv(x, 1) and lv(z, 0))),
             lambda x, z: int(lv(x, 1) and lv(z, 1))]
        P = [lambda r, c: int(lv(r, 2) and lv(c, 0)), lambda r, c: int(lv(r, 2) and lv(c, 1)),
             lambda r, c: int(r == c and lv(c, 0)), lambda r, c: int(lv(r, 0) and lv(c, 0) and r != c),
             lambda r, c: int(lv(r, 0) and lv(c, 1)), lambda r, c: int(lv(r, 1) and lv(c, 0)),
             lambda r, c: int(r == c and lv(c, 1)), lambda r, c: int(lv(r, 1) and lv(c, 1) and r != c)]
        Qv = [lambda x, z: int(lv(x, 0) and lv(z, 0)),
              lambda x, z: int((lv(x, 0) and lv(z, 1)) or (lv(x, 1) and lv(z, 0))),
              lambda x, z: int(lv(x, 1) and lv(z, 1))]
        ker = [[1, 0, 1, 1, 0, 1, 0, 0], [0, 1, 0, 0, 1, 0, 1, 1]]
        names = (['R_t0', 'R_t1', 'R_00', 'R_01', 'R_11'],
                 ['P_t0', 'P_t1', 'P_00d', 'P_00', 'P_01', 'P_10', 'P_11d', 'P_11'], ['Q_00', 'Q_01', 'Q_11'])
        return named, R, P, Qv, ker, [2, 6], names
    if part == 'std':
        named = [T_LAB, nm('a', 0), nm('b', 0)]
        f = lambda x: fo(x, 'a', 'b')
        l0 = lambda x: x[0] == 0
        l1 = lambda x: x[0] == 1
        R = [lambda x, z: (f(z) if x[0] == 2 else 0) + (f(x) if z[0] == 2 else 0),
             lambda x, z: (f(x) + f(z)) * int(l0(x) and l0(z)),
             lambda x, z: f(x) * int(l1(z)) + f(z) * int(l1(x))]
        P = [lambda r, c: f(c) * int(r[0] == 2), lambda r, c: f(r) * int(r == c),
             lambda r, c: f(r) * int(l0(r) and l0(c) and r != c), lambda r, c: f(c) * int(l0(r) and l0(c) and r != c),
             lambda r, c: f(r) * int(l1(c)), lambda r, c: f(c) * int(l1(r))]
        Qv = [lambda x, z: (f(x) + f(z)) * int(l0(x) and l0(z)),
              lambda x, z: f(x) * int(l1(z)) + f(z) * int(l1(x))]
        ker = [[1, 1, 0, 1, 0, 1]]
        names = (['R_t0f', 'R_00f', 'R_01f'], ['P_t0f', 'P_ddf', 'P_00fr', 'P_00fc', 'P_01fr', 'P_10fc'],
                 ['Q_00f', 'Q_01f'])
        return named, R, P, Qv, ker, [1], names
    if part == 'two':
        named = [T_LAB] + [nm(n, 0) for n in 'abcd']
        hv = {frozenset('ac'): 1, frozenset('bd'): 1, frozenset('ad'): -1, frozenset('bc'): -1}

        def h(x, y):
            if x[1] is None or y[1] is None or x == y or x[0] != 0 or y[0] != 0:
                return 0
            return hv.get(frozenset((x[1], y[1])), 0)
        return named, [lambda x, z: h(x, z)], [lambda r, c: h(r, c)], [lambda x, z: h(x, z)], [], [], \
            (['R_h'], ['P_h'], ['Q_h'])
    if part == 'anti':
        named = [T_LAB] + [nm(n, 0) for n in 'abc']
        f1 = lambda x: fo(x, 'a', 'b')
        f2 = lambda x: fo(x, 'b', 'c')
        return named, [], [lambda r, c: f1(r) * f2(c) - f2(r) * f1(c)], [], [], [], ([], ['P_A'], [])
    if part == 'mix':
        named = [T_LAB, nm('a', 0), nm('b', 0), nm('c', 1), nm('d', 1)]
        f = lambda x: fo(x, 'a', 'b')
        g = lambda x: fo(x, 'c', 'd')
        R = [lambda x, z: f(x) * g(z) + f(z) * g(x)]
        P = [lambda r, c: f(r) * g(c), lambda r, c: g(r) * f(c)]
        Qv = [lambda x, z: f(x) * g(z) + f(z) * g(x)]
        return named, R, P, Qv, [], [], (['R_fg'], ['P_fg', 'P_gf'], ['Q_fg'])
    raise ValueError(part)


def keep_for2(part, u0, u1):
    """indices of the R-, P- and Q-vectors forming a basis at level sizes u0, u1 (ints; 3 stands for every size >= 3),
    the gauge P-indices within the kept P-list, and the kernel vectors in the kept basis."""
    if part == 'inv':
        kR = [i for i, need in enumerate([u0 >= 1, u1 >= 1, u0 >= 2, u0 >= 1 and u1 >= 1, u1 >= 2]) if need]
        kP = [i for i, need in enumerate([u0 >= 1, u1 >= 1, u0 >= 1, u0 >= 2, u0 >= 1 and u1 >= 1, u0 >= 1 and u1 >= 1,
                                          u1 >= 1, u1 >= 2]) if need]
        kQ = [i for i, need in enumerate([u0 >= 2, u0 >= 1 and u1 >= 1, u1 >= 2]) if need]
        full = [[1, 0, 1, 1, 0, 1, 0, 0], [0, 1, 0, 0, 1, 0, 1, 1]]
        ker = [[kv[i] for i in kP] for kv in full]
        ker = [v for v in ker if any(v)]
        gauge = [kP.index(i) for i in (2, 6) if i in kP]
        return kR, kP, kQ, gauge, ker
    if part == 'std':      # needs u0 >= 2; at u0 = 2: P3 = -P2 (indices 2, 3), R1 = Q0 = 0
        kR = [0] + ([1] if u0 >= 3 else []) + ([2] if u1 >= 1 else [])
        kP = [0, 1, 2] + ([3] if u0 >= 3 else []) + ([4, 5] if u1 >= 1 else [])
        kQ = ([0] if u0 >= 3 else []) + ([1] if u1 >= 1 else [])
        full = {0: 1, 1: 1, 2: 0, 3: 1, 4: 0, 5: 1}
        v = [full[i] for i in kP]
        if 3 not in kP:
            v[kP.index(2)] = -1
        return kR, kP, kQ, [kP.index(1)], [v]
    if part == 'two':
        return [0], [0], [0], [], []
    if part == 'anti':
        return [], [0], [], [], []
    if part == 'mix':
        return [0], [0, 1], [0], [], []
    raise ValueError(part)


# ---------------------------------------------------------------- pattern sums (levels 0, 1 anonymous; t named)
def patterns2(nslots, named, distinct_pairs):
    out = []

    def rec(i, cur, nan):
        if i == nslots:
            if all(cur[p] != cur[q] for p, q in distinct_pairs):
                out.append(tuple(cur))
            return
        opts = list(named) + [x for x in dict.fromkeys(cur) if x[1] is None]
        for x in opts:
            cur.append(x)
            rec(i + 1, cur, nan)
            cur.pop()
        for lev in (0, 1):
            cur.append((lev, None, ('a', nan)))
            rec(i + 1, cur, nan + 1)
            cur.pop()
    rec(0, [], 0)
    return out


_PAT = {}


def pats(nslots, named, dp):
    key = (nslots, tuple(named), tuple(dp))
    if key not in _PAT:
        _PAT[key] = patterns2(nslots, named, dp)
    return _PAT[key]


def pair_sums2(M, named, vx, vy, kx, ky, fun):
    """R[i][j] = sum over states x (kind kx) and y (kind ky) of vx[i](x) fun(x, y) vy[j](y), by patterns.
    kinds: 'pair' (unordered distinct pair, enumerated ordered with weight 1/2) or 'ord' (ordered, may coincide)."""
    dp = []
    if kx == 'pair':
        dp.append((0, 1))
    if ky == 'pair':
        dp.append((2, 3))
    w = Fr(1, 2 ** len(dp))
    nnamed = [sum(1 for x in named if x[0] == i) for i in (0, 1)]
    acc = {}
    for labs in pats(4, named, dp):
        cx = [f(labs[0], labs[1]) for f in vx]
        if not any(cx):
            continue
        cy = [g(labs[2], labs[3]) for g in vy]
        if not any(cy):
            continue
        anon = set(x for x in labs if x[1] is None)
        mv = tuple(sum(1 for x in anon if x[0] == i) for i in (0, 1))
        key = (canon(labs), mv)
        A = acc.setdefault(key, [[0] * len(vy) for _ in vx])
        for i, a in enumerate(cx):
            if a:
                for j, bb in enumerate(cy):
                    if bb:
                        A[i][j] += a * bb
    one = M.one
    byval = {}
    for (vk, mv), A in acc.items():
        cnt = ff(M.U[0] - nnamed[0], mv[0], one) * ff(M.U[1] - nnamed[1], mv[1], one)
        if cnt == 0:
            continue
        B = byval.setdefault(vk, [[M.zero] * len(vy) for _ in vx])
        for i in range(len(vx)):
            for j in range(len(vy)):
                if A[i][j]:
                    B[i][j] = B[i][j] + cnt * A[i][j]
    R = [[M.zero] * len(vy) for _ in vx]
    for vk, B in byval.items():
        val = fun(*from_key(vk))
        if val == 0:
            continue
        for i in range(len(vx)):
            for j in range(len(vy)):
                if B[i][j] != 0:
                    R[i][j] = R[i][j] + val * B[i][j]
    wf = num(one, w)
    return [[wf * v for v in row] for row in R]


def single_sums2(M, named, vx, vy, kind, valfun):
    """G[i][j] = sum over states x of kind of vx[i](x) vy[j](x) valfun(x)."""
    dp = [(0, 1)] if kind == 'pair' else []
    w = Fr(1, 2 ** len(dp))
    nnamed = [sum(1 for x in named if x[0] == i) for i in (0, 1)]
    one = M.one
    G = [[M.zero] * len(vy) for _ in vx]
    for labs in pats(2, named, dp):
        cx = [f(labs[0], labs[1]) for f in vx]
        cy = [f(labs[0], labs[1]) for f in vy]
        if not any(cx) or not any(cy):
            continue
        anon = set(x for x in labs if x[1] is None)
        mv = [sum(1 for x in anon if x[0] == i) for i in (0, 1)]
        cnt = ff(M.U[0] - nnamed[0], mv[0], one) * ff(M.U[1] - nnamed[1], mv[1], one)
        if cnt == 0:
            continue
        v = valfun(*labs) * cnt
        for i in range(len(vx)):
            for j in range(len(vy)):
                if cx[i] * cy[j]:
                    G[i][j] = G[i][j] + v * (cx[i] * cy[j])
    wf = num(one, w)
    return [[wf * v for v in row] for row in G]


def reduce_part2(M, part, u0, u1):
    """reduced forms of the part at the model M (sizes u0, u1 as ints for the basis choice: 3 = every size >= 3).
    Returns dict with the blocks RR, RP, RQ, PP, PQ, QQ (for the model's bscale), the Gram matrices, names, kernel
    vectors (on the P-vectors) and the gauge indices."""
    named, R, P, Qv, ker0, gauge0, names = part_spec2(part)
    kR, kP, kQ, gauge, ker = keep_for2(part, u0, u1)
    R = [R[i] for i in kR]
    P = [P[i] for i in kP]
    Qv = [Qv[i] for i in kQ]
    names = ([names[0][i] for i in kR], [names[1][i] for i in kP], [names[2][i] for i in kQ])
    one, zero = M.one, M.zero
    s = M.s
    out = dict(names=names, ker=ker, gauge=gauge)
    out['GP'] = single_sums2(M, named, P, P, 'ord', lambda r, c: one)
    out['PP'] = pair_sums2(M, named, P, P, 'ord', 'ord', M.PP)
    if R:
        out['GR'] = single_sums2(M, named, R, R, 'pair', lambda x, z: one)
        out['RR'] = single_sums2(M, named, R, R, 'pair', lambda x, z: M.RR((x, z), (x, z)))
        out['RP'] = pair_sums2(M, named, R, P, 'pair', 'ord', M.RP)
    if Qv:
        out['GQ'] = single_sums2(M, named, Qv, Qv, 'pair', lambda x, z: one)
        out['DQ'] = single_sums2(M, named, Qv, Qv, 'pair', lambda x, z: s(x) * s(z))
        out['PQ'] = pair_sums2(M, named, P, Qv, 'ord', 'pair', M.PQ)
        RW = pair_sums2(M, named, Qv, Qv, 'pair', 'pair', M.W4)
        GQ, DQ = out['GQ'], out['DQ']
        full = mm(mm(GQ, minv(RW, one), zero), GQ, zero)
        out['QQ'] = [[2 * DQ[i][j] + M.b * (2 * full[i][j] - 2 * DQ[i][j]) for j in range(len(Qv))]
                     for i in range(len(Qv))]
        if R:
            # RQ = -E~ + 1/2 G_R RA2^-1 RN G_Q^-1 DQ   (E~: the embedding of the pairs of U, rescaled: s_P)
            RE = single_sums2(M, named, R, Qv, 'pair', lambda x, z: s(x) * s(z))
            RA2 = pair_sums2(M, named, R, R, 'pair', 'pair', M.A2hat)
            RN = pair_sums2(M, named, R, Qv, 'pair', 'pair', M.N2c)
            T = mm(mm(mm(mm(out['GR'], minv(RA2, one), zero), RN, zero), minv(GQ, one), zero), DQ, zero)
            out['RQ'] = [[-RE[i][j] + M.b * (T[i][j] / 2) for j in range(len(Qv))] for i in range(len(R))]
    return out


def assemble(red):
    """the reduced form on (R-vectors, P-vectors, Q-vectors) as one symmetric matrix, and the index lists."""
    nR = len(red['names'][0])
    nP = len(red['names'][1])
    nQ = len(red['names'][2])
    n = nR + nP + nQ
    zero = red['PP'][0][0] * 0
    F = [[zero] * n for _ in range(n)]

    def put(blk, i0, j0, ni, nj):
        for i in range(ni):
            for j in range(nj):
                F[i0 + i][j0 + j] = blk[i][j]
                F[j0 + j][i0 + i] = blk[i][j]
    if nR:
        put(red['RR'], 0, 0, nR, nR)
        put(red['RP'], 0, nR, nR, nP)
    put(red['PP'], nR, nR, nP, nP)
    if nQ:
        put(red['PQ'], nR, nR + nP, nP, nQ)
        put(red['QQ'], nR + nP, nR + nP, nQ, nQ)
        if nR:
            put(red['RQ'], 0, nR + nP, nR, nQ)
    return F, list(range(nR)), list(range(nR, nR + nP)), list(range(nR + nP, n))
