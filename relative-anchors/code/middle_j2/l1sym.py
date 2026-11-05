"""Symmetry reduction of the letter-1 block of (R_inf)_1 at two weight levels, exact, symbolic in the level sizes.

Setting: relative-anchors/README.md (Sections 1, 2) and middle_j1/NOTES.md (Section 4).  The copies form two levels:
u_0 copies of squared weight Q_0 and u_1 copies of squared weight Q_1, k = u_0 + u_1, L = |s| = u_0 Q_0 + u_1 Q_1.
States: P_{r,c} = (1@r, 1@c) and Q_{ab} = (vac, 1@{a,b}).  Units 2^{k-1}, rescaled frame of ylib (the entry between
two states is multiplied by the product of the weights w_t over the letters 1 of both states), in which every entry
is rational in the squared weights.  The closed forms (NOTES.md, Section 2), with s_c = w_c^2, d_c = 2(L - s_c),
g0 = sum_c s_c/d_c, sg = 2 + g0:

  PP: Y~[(r1,c1),(r2,c2)] = (1/2 - [r1=r2])[c2=r1][c1=r2] s_r1 s_r2 + [r1=r2][c1=c2](2 - [r1=c1]) s_r1 s_c1
                            - b * beta~((r2,c1),(r1,c2)),
      beta~((r,c),(r',c')) = s_r s_r' g0/(2 sg)                        if c = r, c' = r',
                           = 2 s_r s_r' s_c'/(sg d_c')                  if c = r, c' != r'  (and symmetrically),
                           = 4 s_r s_r' (s_c [c=c']/d_c - s_c s_c'/(d_c d_c' sg))   otherwise;
  PQ: Y~[(r,c),{a,b}] = -s_r Ach[c, c*]                                if r is in {a,b}, c* the other copy,
                      = s_r (s_a Ach[c,b] + s_b Ach[c,a])/(L - s_a - s_b)   otherwise,
      with Ach = Ahat^{-1}, Ahat = D^{-1}(I + N1/2)D^{-1}: Ahat[c,c] = (1 + N1[c,c]/2)/s_c,
      Ahat[c,c'] = 1/(4(L - s_c - s_c')), N1[c,c] = sum_{u != c} s_u/(2(L - s_c - s_u))
      (for bscale b: A-part -[r in P][c = c*] s_a s_b plus b times the rest);
  QQ: Y~_QQ = 2 W^{-1} with W[P,P] = (1 + N2[P,P]/2)/(s_a s_b), N2[P,P] = sum_{u not in P} s_u/(2(L - s_a - s_b - s_u)),
      W[{c,a},{c,b}] = 1/(4 s_c (L - s_a - s_b - s_c)) (a != b), 0 for disjoint pairs
      (for bscale b: Y~_QQ = 2 D2 + b (2 W^{-1} - 2 D2), D2 = diag(s_a s_b)).

Symmetry: S_{u_0} x S_{u_1} permutes the copies inside each level.  For each irrep type a short list of explicit test
vectors (equivariant images of one seed vector) spans the multiplicity space; the reduced forms are
R[i][j] = <v_i, Y~ v_j>, computed by summing over the joint patterns of the indices of two states (each index a named
copy of the seed, or an anonymous copy of level 0 or 1; anonymous copies are distinct from the named ones and from each
other unless equal), weighted by the number of realisations (falling factorials in u_i).  Values depend only on the
levels and the equalities of the indices.  Everything is generic over the number type: Fractions for exact checks,
flintfield elements for the symbolic certificates.

Parts (irrep on level 0 unless stated; f = e_a - e_b with named copies a, b of level 0):
  inv   (trivial):       P: T1 = 1[r,c in U0, r != c], T2 = 1[r = c in U0], T3 = 1[r in U0, c in U1], T4 = 1[r in U1, c in U0],
                         T5 = 1[r,c in U1, r != c], T6 = 1[r = c in U1];  Q: T7 = 1[a,b in U0], T8 = 1[across], T9 = 1[a,b in U1].
                         kernel: T1 + T2 + T4 and T3 + T5 + T6.
  std   (standard of level 0):  P: S1 = f_r [r,c in U0, r != c], S2 = f_c [same], S3 = f_r [r = c], S4 = f_r [c in U1],
                         S5 = f_c [r in U1];  Q: S6 = f_a + f_b [a,b in U0], S7 = f_a [b in U1] + f_b [a in U1].
                         kernel: S2 + S3 + S5.
  two   ((u0-2,2) of level 0): h = e_ac + e_bd - e_ad - e_bc on pairs of four named copies a,b,c,d of level 0;
                         P: H1 = h_{rc} (r != c);  Q: H2 = h.
  anti  ((u0-2,1,1) of level 0): A = f1 ^ f2, f1 = e_a - e_b, f2 = e_b - e_c;  P: A1 = A[r,c].
  mix   (standard of level 0 times standard of level 1): f = e_a - e_b (level 0), g = e_c - e_d (level 1);
                         P: B1 = f_r g_c, B2 = g_r f_c;  Q: B3 = f_a g_b + f_b g_a.
"""
from fractions import Fraction as Fr


# ---------------------------------------------------------------- the model (closed forms)
def num(one, v):
    """the rational number v in the number type of one (Fraction or a flintfield element)."""
    v = Fr(v)
    if isinstance(one, Fr):
        return v
    import flint
    return one.K(flint.fmpq(v.numerator, v.denominator))


class Model:
    def __init__(self, one, Q, U, b=1):
        """one: the unit of the number type; Q = [Q0, Q1] squared weights; U = [u0, u1] level sizes (numbers of the
        same type); b: B-scale (1 for the true block)."""
        self.one = one
        self.zero = one * 0
        self.Q, self.U = list(Q), list(U)
        self.b = num(one, b)
        L = U[0] * Q[0] + U[1] * Q[1]
        self.L = L
        self.g0 = U[0] * Q[0] / (2 * (L - Q[0])) + U[1] * Q[1] / (2 * (L - Q[1]))
        self.sg = self.g0 + 2
        # Ahat and its inverse Ach (invariant k x k matrices: diagonal d[i], off-diagonal e[i][j] by levels)
        N1 = [sum(((U[j] - (1 if i == j else 0)) * Q[j] / (2 * (L - Q[i] - Q[j])) for j in (0, 1)), self.zero)
              for i in (0, 1)]
        d = [(one + N1[i] / 2) / Q[i] for i in (0, 1)]
        e = [[one / (4 * (L - Q[i] - Q[j])) for j in (0, 1)] for i in (0, 1)]
        lam = [one / (d[i] - e[i][i]) for i in (0, 1)]
        Mi = [[d[0] + (U[0] - 1) * e[0][0], U[1] * e[0][1]], [U[0] * e[1][0], d[1] + (U[1] - 1) * e[1][1]]]
        dt = Mi[0][0] * Mi[1][1] - Mi[0][1] * Mi[1][0]
        m = [[Mi[1][1] / dt, -Mi[0][1] / dt], [-Mi[1][0] / dt, Mi[0][0] / dt]]
        # (a level of size 0 is allowed; its entries are never used, since no copy of that level exists)
        dv = lambda x, u: self.zero if u == 0 else x / u
        ep = [[dv(m[0][0] - lam[0], U[0]), dv(m[0][1], U[1])], [dv(m[1][0], U[0]), dv(m[1][1] - lam[1], U[1])]]
        self.Ad = [lam[i] + ep[i][i] for i in (0, 1)]
        self.Ae = ep
        self.Ahat = (d, e)
        self._cache = {}

    def s(self, x):
        return self.Q[x[0]]

    def dl(self, x):
        return 2 * (self.L - self.s(x))

    def Ach(self, c, cp):
        return self.Ad[c[0]] if c == cp else self.Ae[c[0]][cp[0]]

    def beta(self, r, c, rp, cp):
        s = self.s
        if c == r and cp == rp:
            return s(r) * s(rp) * self.g0 / (2 * self.sg)
        if c == r:
            return 2 * s(r) * s(rp) * s(cp) / (self.sg * self.dl(cp))
        if cp == rp:
            return 2 * s(rp) * s(r) * s(c) / (self.sg * self.dl(c))
        out = -s(c) * s(cp) / (self.dl(c) * self.dl(cp) * self.sg)
        if c == cp:
            out = out + s(c) / self.dl(c)
        return 4 * s(r) * s(rp) * out

    def PP(self, r1, c1, r2, c2):
        s = self.s
        A = self.zero
        if c2 == r1 and c1 == r2:
            A = A + (self.one / 2 - (1 if r1 == r2 else 0)) * s(r1) * s(r2)
        if r1 == r2 and c1 == c2:
            A = A + (2 - (1 if r1 == c1 else 0)) * s(r1) * s(c1)
        return A - self.b * self.beta(r2, c1, r1, c2)

    def PQ(self, r, c, a, bb):
        s = self.s
        if r == a or r == bb:
            cs = bb if r == a else a
            full = -s(r) * self.Ach(c, cs)
            Ap = -s(a) * s(bb) if c == cs else self.zero
        else:
            full = s(r) * (s(a) * self.Ach(c, bb) + s(bb) * self.Ach(c, a)) / (self.L - s(a) - s(bb))
            Ap = self.zero
        return Ap + self.b * (full - Ap)

    def W(self, a, bb, a2, b2):
        s, U, Q, L = self.s, self.U, self.Q, self.L
        P, P2 = {a, bb}, {a2, b2}
        if P == P2:
            n2 = self.zero
            for i in (0, 1):
                cnt = sum(1 for x in P if x[0] == i)
                n2 = n2 + (U[i] - cnt) * Q[i] / (2 * (L - s(a) - s(bb) - Q[i]))
            return (self.one + n2 / 2) / (s(a) * s(bb))
        I = P & P2
        if len(I) == 1:
            c = next(iter(I))
            x = next(iter(P - I))
            y = next(iter(P2 - I))
            return self.one / (4 * s(c) * (L - s(x) - s(y) - s(c)))
        return self.zero

    def D2(self, a, bb, a2, b2):
        return self.s(a) * self.s(bb) if {a, bb} == {a2, b2} else self.zero


# ---------------------------------------------------------------- test vectors
def lab_named(name, lev):
    return (lev, name, name)


def f_of(x, plus, minus):
    return 1 if x[1] == plus else (-1 if x[1] == minus else 0)


def part_spec(part):
    """returns (named labels, P-vectors, Q-vectors, kernel vectors (coefficient lists on the P-vectors), names).
    The irrep sits on level 0; for 'mix' on both levels."""
    if part == 'inv':
        named = []
        P = [lambda r, c: int(r[0] == 0 and c[0] == 0 and r != c), lambda r, c: int(r == c and r[0] == 0),
             lambda r, c: int(r[0] == 0 and c[0] == 1), lambda r, c: int(r[0] == 1 and c[0] == 0),
             lambda r, c: int(r[0] == 1 and c[0] == 1 and r != c), lambda r, c: int(r == c and r[0] == 1)]
        Qv = [lambda a, b: int(a[0] == 0 and b[0] == 0), lambda a, b: int(a[0] != b[0]),
              lambda a, b: int(a[0] == 1 and b[0] == 1)]
        ker = [[1, 1, 0, 1, 0, 0], [0, 0, 1, 0, 1, 1]]
        names = ['T1', 'T2', 'T3', 'T4', 'T5', 'T6'], ['T7', 'T8', 'T9']
        return named, P, Qv, ker, names
    if part == 'std':
        named = [lab_named('a', 0), lab_named('b', 0)]
        f = lambda x: f_of(x, 'a', 'b')
        P = [lambda r, c: f(r) * int(r[0] == 0 and c[0] == 0 and r != c),
             lambda r, c: f(c) * int(r[0] == 0 and c[0] == 0 and r != c),
             lambda r, c: f(r) * int(r == c),
             lambda r, c: f(r) * int(c[0] == 1),
             lambda r, c: f(c) * int(r[0] == 1)]
        Qv = [lambda a, b: (f(a) + f(b)) * int(a[0] == 0 and b[0] == 0),
              lambda a, b: f(a) * int(b[0] == 1) + f(b) * int(a[0] == 1)]
        ker = [[0, 1, 1, 0, 1]]
        names = ['S1', 'S2', 'S3', 'S4', 'S5'], ['S6', 'S7']
        return named, P, Qv, ker, names
    if part == 'two':
        named = [lab_named(n, 0) for n in 'abcd']
        hv = {frozenset('ac'): 1, frozenset('bd'): 1, frozenset('ad'): -1, frozenset('bc'): -1}

        def h(x, y):
            if x[1] is None or y[1] is None or x == y:
                return 0
            return hv.get(frozenset((x[1], y[1])), 0)
        P = [lambda r, c: h(r, c)]
        Qv = [lambda a, b: h(a, b)]
        return named, P, Qv, [], (['H1'], ['H2'])
    if part == 'anti':
        named = [lab_named(n, 0) for n in 'abc']
        f1 = lambda x: f_of(x, 'a', 'b')
        f2 = lambda x: f_of(x, 'b', 'c')
        P = [lambda r, c: f1(r) * f2(c) - f2(r) * f1(c)]
        return named, P, [], [], (['A1'], [])
    if part == 'mix':
        named = [lab_named('a', 0), lab_named('b', 0), lab_named('c', 1), lab_named('d', 1)]
        f = lambda x: f_of(x, 'a', 'b')
        g = lambda x: f_of(x, 'c', 'd')
        P = [lambda r, c: f(r) * g(c), lambda r, c: g(r) * f(c)]
        Qv = [lambda a, b: f(a) * g(b) + f(b) * g(a)]
        return named, P, Qv, [], (['B1', 'B2'], ['B3'])
    raise ValueError(part)


# ---------------------------------------------------------------- pattern sums
def patterns(nslots, named, distinct_pairs):
    """all assignments of nslots index slots to labels: a named label, or an anonymous label (lev, None, ('a', i));
    anonymous labels are numbered in order of first use.  distinct_pairs: list of slot pairs that must differ."""
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


def ff(n, m, one):
    out = one
    for i in range(m):
        out = out * (n - i)
    return out


def canon(labs):
    idx = {}
    key = []
    for x in labs:
        if x not in idx:
            idx[x] = len(idx)
        key.append((x[0], idx[x]))
    return tuple(key)


def from_key(key):
    return [(lev, None, i) for lev, i in key]


def pair_sums(M, named, vx, vy, kx, ky, valname):
    """R[i][j] = sum over states x of kind kx and y of kind ky of vx[i](x) val(x, y) vy[j](y), summed by patterns.
    Q-states are unordered pairs: enumerated as ordered pairs with weight 1/2 each."""
    dp = []
    if kx == 'Q':
        dp.append((0, 1))
    if ky == 'Q':
        dp.append((2, 3))
    w = Fr(1, 2 ** len(dp))
    pats = patterns(4, named, dp)
    nnamed = [sum(1 for x in named if x[0] == i) for i in (0, 1)]
    acc = {}
    for labs in pats:
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
    fun = getattr(M, valname)
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
    wf = one * w.numerator / w.denominator
    for vk, B in byval.items():
        val = fun(*from_key(vk))
        if val == 0:
            continue
        for i in range(len(vx)):
            for j in range(len(vy)):
                if B[i][j] != 0:
                    R[i][j] = R[i][j] + val * B[i][j]
    return [[wf * v for v in row] for row in R]


def single_sums(M, named, vx, kind, valfun):
    """G[i][j] = sum over states x of kind of vx[i](x) vx[j](x) valfun(x) (valfun on labels)."""
    dp = [(0, 1)] if kind == 'Q' else []
    w = Fr(1, 2 ** len(dp))
    pats = patterns(2, named, dp)
    nnamed = [sum(1 for x in named if x[0] == i) for i in (0, 1)]
    one = M.one
    G = [[M.zero] * len(vx) for _ in vx]
    for labs in pats:
        cx = [f(labs[0], labs[1]) for f in vx]
        if not any(cx):
            continue
        anon = set(x for x in labs if x[1] is None)
        mv = [sum(1 for x in anon if x[0] == i) for i in (0, 1)]
        cnt = ff(M.U[0] - nnamed[0], mv[0], one) * ff(M.U[1] - nnamed[1], mv[1], one)
        if cnt == 0:
            continue
        v = valfun(*labs) * cnt
        for i in range(len(vx)):
            for j in range(len(vx)):
                if cx[i] * cx[j]:
                    G[i][j] = G[i][j] + v * (cx[i] * cx[j])
    wf = one * w.numerator / w.denominator
    return [[wf * v for v in row] for row in G]


# ---------------------------------------------------------------- small matrix algebra
def mm(A, B, zero):
    return [[sum((A[i][m] * B[m][j] for m in range(len(B))), zero) for j in range(len(B[0]))] for i in range(len(A))]


def mT(A):
    return [list(r) for r in zip(*A)]


def minv(A, one):
    n = len(A)
    zero = one * 0
    M = [list(A[i]) + [one if i == j else zero for j in range(n)] for i in range(n)]
    for c in range(n):
        r = next(i for i in range(c, n) if M[i][c] != 0)
        M[c], M[r] = M[r], M[c]
        pv = M[c][c]
        M[c] = [v / pv for v in M[c]]
        for i in range(n):
            if i != c and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * bb for a, bb in zip(M[i], M[c])]
    return [row[n:] for row in M]


def det(A, one):
    A = [list(r) for r in A]
    n = len(A)
    d = one
    for c in range(n):
        r = next((i for i in range(c, n) if A[i][c] != 0), None)
        if r is None:
            return one * 0
        if r != c:
            A[c], A[r] = A[r], A[c]
            d = -d
        pv = A[c][c]
        d = d * pv
        for i in range(c + 1, n):
            if A[i][c] != 0:
                f = A[i][c] / pv
                A[i] = [a - f * bb for a, bb in zip(A[i], A[c])]
    return d


# ---------------------------------------------------------------- the reduced forms of one part
def reduce_part(M, part, keepP=None, keepQ=None):
    """reduced forms of the part at the model M.  keepP/keepQ: indices of the test vectors to keep (drop vectors that
    vanish or are dependent at small level sizes).  Returns dict with RPP, RPQ, RW, GQ, RD, RQQ (the reduced vacuum
    block, for the model's bscale), S (Schur complement onto the P-vectors) and the kernel vectors."""
    named, P, Qv, ker, names = part_spec(part)
    if keepP is not None:
        P = [P[i] for i in keepP]
        ker = kernel_coeffs(part, keepP)
        names = ([names[0][i] for i in keepP], names[1])
    if keepQ is not None:
        Qv = [Qv[i] for i in keepQ]
        names = (names[0], [names[1][i] for i in keepQ])
    one, zero = M.one, M.zero
    out = dict(names=names, ker=ker)
    out['RPP'] = pair_sums(M, named, P, P, 'P', 'P', 'PP')
    out['GP'] = single_sums(M, named, P, 'P', lambda r, c: one)
    if Qv:
        out['RPQ'] = pair_sums(M, named, P, Qv, 'P', 'Q', 'PQ')
        out['RW'] = pair_sums(M, named, Qv, Qv, 'Q', 'Q', 'W')
        out['GQ'] = single_sums(M, named, Qv, 'Q', lambda a, b: one)
        out['RD'] = single_sums(M, named, Qv, 'Q', lambda a, b: M.s(a) * M.s(b))
        GQ, RW, RD = out['GQ'], out['RW'], out['RD']
        RWi = minv(RW, one)
        Rinv = mm(mm(GQ, RWi, zero), GQ, zero)                 # reduced form of W^{-1}
        out['RQQ'] = [[2 * RD[i][j] + M.b * (2 * Rinv[i][j] - 2 * RD[i][j]) for j in range(len(Qv))]
                      for i in range(len(Qv))]
        RPQ = out['RPQ']
        if M.b == 1:
            GQi = minv(GQ, one)
            T = mm(mm(mm(RPQ, GQi, zero), RW, zero), GQi, zero)
            C = mm(T, mT(RPQ), zero)
            out['S'] = [[out['RPP'][i][j] - C[i][j] / 2 for j in range(len(P))] for i in range(len(P))]
        else:
            Qi = minv(out['RQQ'], one)
            C = mm(mm(RPQ, Qi, zero), mT(RPQ), zero)
            out['S'] = [[out['RPP'][i][j] - C[i][j] for j in range(len(P))] for i in range(len(P))]
    else:
        out['S'] = out['RPP']
    return out


# regimes: which test vectors exist at small level sizes (u_i given as an int, or None when symbolic >= 3)
def keep_for(part, u0, u1):
    """indices of the P- and Q-vectors that form a basis at level sizes u0, u1 (ints; a value >= 3 stands for every
    size >= 3; u1 = 0 is the one-level case).  Returns (keepP, keepQ, gauge), gauge = the P-indices (within keepP) of
    the vectors supported on the diagonal states P_{c,c}, which are deleted to fix the kernel."""
    if part == 'inv':
        kp = [i for i, need in enumerate([u0 >= 2, u0 >= 1, u0 >= 1 and u1 >= 1, u0 >= 1 and u1 >= 1, u1 >= 2,
                                          u1 >= 1]) if need]
        kq = [i for i, need in enumerate([u0 >= 2, u0 >= 1 and u1 >= 1, u1 >= 2]) if need]
        return kp, kq, [kp.index(i) for i in (1, 5) if i in kp]
    if part == 'std':        # at u0 = 2: S2 = -S1 and S6 = 0
        kp = [0] + ([1] if u0 >= 3 else []) + [2] + ([3, 4] if u1 >= 1 else [])
        kq = ([0] if u0 >= 3 else []) + ([1] if u1 >= 1 else [])
        return kp, kq, [kp.index(2)]
    if part == 'two':
        return [0], [0], []
    if part == 'anti':
        return [0], None, []
    if part == 'mix':
        return [0, 1], [0], []
    raise ValueError(part)


def kernel_coeffs(part, keepP):
    """kernel vectors of the reduced P-form in the kept basis: the S_{u0} x S_{u1}-parts of the indicators 1_c of the
    column groups (in rescaled coordinates the kernel vectors p (x) |1@c> of the block)."""
    if part == 'inv':
        full = [[1, 1, 0, 1, 0, 0], [0, 0, 1, 0, 1, 1]]
        out = [[kv[i] for i in keepP] for kv in full]
        return [v for v in out if any(v)]
    if part == 'std':
        full = {0: 0, 1: 1, 2: 1, 3: 0, 4: 1}
        v = [full[i] for i in keepP]
        if 1 not in keepP:   # u0 = 2: S2 = -S1
            v[keepP.index(0)] = -1
        return [v]
    return []
