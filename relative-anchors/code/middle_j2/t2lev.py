"""Family 1: the type I blocks of (R_inf)_2 with T = {t, tau} (one letter 1), weights arbitrary on t and tau, and the
U-weights on l levels, symbolically in the level sizes.

Copies: t (squared weight a), tau (squared weight b), and U = U_0 + ... + U_{l-1}, |U_i| = u_i, squared weight q_i on
U_i (q_0 = 1 in the certificates: every rescaled entry is homogeneous of degree 1 in the squared weights).
L = a + b + sum u_i q_i.  Units 2^{k-2}.  States (NOTES.md, Section 1):
  R_c  = (2@t,   2@tau 1@c), c in Ut = U + {t};     R'_c = (2@tau, 2@t 1@c), c in UT = U + {tau};
  P_r  = (1@r,   2@t 2@tau), r in [k];             Q_c  = (vac,   2@t 2@tau 1@c), c in U.
Rescaled frame D = diag(mu), mu = w_c on every state (the copy carrying the letter 1); O~ = D O D.

Invariant-operator algebra: the group S_{u_0} x ... x S_{u_{l-1}} permutes the copies inside each level; t and tau are
fixed.  A single-index space is a tuple of orbits ('t', 'T' for tau, and level numbers).  An invariant operator between
two such spaces is stored by its reduced matrix M[o_y][o_x] = sum_{x in o_x} O[y, x] (y in o_y fixed) and by its scalar
sd[i] = O[h, h] - O[h, h'] on the standard part of each level i (h != h' in U_i).  Products, transposes (Gram matrices
diag(orbit sizes)) and inverses are small computations.  Rescaled products: (O1 O2)~ = O1~ D^-2 O2~,
(O^-1)~ = D^2 (O~)^-1 D^2.

build(...) returns the invariant form and the standard forms of
  Phi = C - B^T diag(1/N) B on (R, R', Q)   (Lemma P2; only for bscale = 1), and of
  the gauge-fixed block G (the state P_t deleted) on (R, R', P without t, Q), with the B-part scaled by bscale.
"""
from fractions import Fraction as Fr
import flint
from flintfield import Field


def make(K, a, b, qs, us, bscale=1):
    one, zero = K(1), K(0)
    l = len(qs)
    LV = list(range(l))

    def kf(v):
        v = Fr(v)
        return K(flint.fmpq(v.numerator, v.denominator))

    bs = kf(bscale)
    size = {'t': one, 'T': one}
    sq = {'t': a, 'T': b}
    for i in LV:
        size[i] = us[i]
        sq[i] = qs[i]
    L = a + b + sum((us[i] * qs[i] for i in LV), zero)

    def mm(A, B):
        return [[sum((A[i][m] * B[m][j] for m in range(len(B))), zero) for j in range(len(B[0]))] for i in range(len(A))]

    def mT(A):
        return [list(r) for r in zip(*A)]

    def minv(A):
        n = len(A)
        M = [list(A[i]) + [one if i == j else zero for j in range(n)] for i in range(n)]
        for c in range(n):
            r = next(i for i in range(c, n) if M[i][c] != 0)
            M[c], M[r] = M[r], M[c]
            pv = M[c][c]
            M[c] = [v / pv for v in M[c]]
            for i in range(n):
                if i != c and M[i][c] != 0:
                    f = M[i][c]
                    M[i] = [x - f * y for x, y in zip(M[i], M[c])]
        return [row[n:] for row in M]

    def G(space):
        return [[size[o] if i == j else zero for j, o in enumerate(space)] for i, o2 in enumerate(space)]

    def Gi(space):
        return [[one / size[o] if i == j else zero for j, o in enumerate(space)] for i, o2 in enumerate(space)]

    class Op:
        def __init__(self, dom, cod, M, sd):
            self.dom, self.cod, self.M, self.sd = tuple(dom), tuple(cod), M, sd

        def __add__(self, o):
            assert self.dom == o.dom and self.cod == o.cod
            return Op(self.dom, self.cod, [[x + y for x, y in zip(r1, r2)] for r1, r2 in zip(self.M, o.M)],
                      [x + y for x, y in zip(self.sd, o.sd)])

        def __sub__(self, o):
            return self + o.scale(-one)

        def scale(self, c):
            c = c if not isinstance(c, (int, Fr)) else kf(c)
            return Op(self.dom, self.cod, [[c * v for v in r] for r in self.M], [c * v for v in self.sd])

        def __matmul__(self, o):
            assert self.dom == o.cod, (self.dom, o.cod)
            return Op(o.dom, self.cod, mm(self.M, o.M), [x * y for x, y in zip(self.sd, o.sd)])

        def T(self):
            return Op(self.cod, self.dom, mm(mm(Gi(self.dom), mT(self.M)), G(self.cod)), list(self.sd))

        def inv(self):
            assert self.dom == self.cod
            return Op(self.dom, self.cod, minv(self.M), [one / v for v in self.sd])

        def restrict(self, dom=None, cod=None):
            dom = tuple(dom) if dom is not None else self.dom
            cod = tuple(cod) if cod is not None else self.cod
            M = [[self.M[self.cod.index(oy)][self.dom.index(ox)] for ox in dom] for oy in cod]
            return Op(dom, cod, M, list(self.sd))

    def op(dom, cod, f):
        """f(oy, ox, same) = the rescaled entry O~[y, x] for y in orbit oy of cod, x in orbit ox of dom; same: y = x."""
        M = []
        for oy in cod:
            row = []
            for ox in dom:
                if oy == ox and oy in ('t', 'T'):
                    row.append(f(oy, ox, True))
                elif oy == ox:
                    # a level of fixed size 1 has no pair of distinct copies: its off-diagonal value is not used
                    off = zero if size[oy] - 1 == 0 else (size[oy] - 1) * f(oy, oy, False)
                    row.append(f(oy, oy, True) + off)
                else:
                    row.append(size[ox] * f(oy, ox, False))
            M.append(row)
        sd = [f(i, i, True) - (f(i, i, False) if size[i] - 1 != 0 else zero) for i in LV]
        return Op(dom, cod, M, sd)

    Kk = ('t', 'T') + tuple(LV)
    Ut = ('t',) + tuple(LV)
    UT = ('T',) + tuple(LV)
    U = tuple(LV)
    s = lambda o: sq[o]
    D2 = {sp: op(sp, sp, lambda oy, ox, same: s(oy) if same else zero) for sp in (Kk, Ut, UT, U)}
    Dm2 = {sp: op(sp, sp, lambda oy, ox, same: one / s(oy) if same else zero) for sp in (Kk, Ut, UT, U)}

    def comp(*ops):
        out = ops[0]
        for o in ops[1:]:
            out = out @ Dm2[out.dom] @ o
        return out

    def rinv(o):
        return D2[o.cod] @ o.inv() @ D2[o.dom]

    def osum(space, excl, fun):
        """sum over the copies x of the orbits of space, x != the copy excl (given by its orbit), of fun(orbit of x)."""
        tot = zero
        for o in space:
            n = size[o] - (1 if o == excl else 0)
            if n != 0:
                tot = tot + n * fun(o)
        return tot

    def eps(flip, o):
        return -one if o == flip else one

    # ---- N_{XY} = J^-_X C^{-1} J^+_Y on a single-index sector, C(x, z) = cfun(x, z) on the pairs ----
    def Nop(space, cfun, X, Y):
        def f(oy, ox, same):
            if same:
                return s(oy) * osum(space, oy, lambda o: eps(X, o) * eps(Y, o) * s(o) / cfun(oy, o))
            return s(oy) * s(ox) * eps(Y, oy) * eps(X, ox) / cfun(oy, ox)
        return op(space, space, f)

    # ---- R-R' ----
    C2 = lambda x, z: 2 * (L - s(x) - s(z))
    A2 = D2[Kk] + Nop(Kk, C2, None, None).scale(Fr(1, 2))
    W = Nop(Kk, C2, 'T', 't') - comp(Nop(Kk, C2, 'T', None), rinv(A2), Nop(Kk, C2, None, 't')).scale(Fr(1, 2))
    W = W.restrict(dom=UT, cod=Ut)
    A_RRp = op(UT, Ut, lambda oy, ox, same: s(oy) / 2 if (same and oy not in ('t', 'T')) else zero)
    Y_RRp = A_RRp - W.scale(Fr(1, 4))
    B_RRp = W.scale(Fr(1, 4))

    # ---- diagonal blocks of R, R' ----
    D_R = op(Ut, Ut, lambda oy, ox, same: s(oy) * (1 if oy == 't' else 2) if same else zero)
    D_Rp = op(UT, UT, lambda oy, ox, same: s(oy) * (1 if oy == 'T' else 2) if same else zero)

    # ---- P-P ----
    Nr = lambda o: 1 if o in ('t', 'T') else 2
    LP = 2 * L - a - b
    A_PP = op(Kk, Kk, lambda oy, ox, same: s(oy) * Nr(oy) if same else zero)
    B_PP = op(Kk, Kk, lambda oy, ox, same: s(oy) * Nr(oy) * s(ox) * Nr(ox) / LP)
    Ninv = op(Kk, Kk, lambda oy, ox, same: s(oy) / Nr(oy) if same else zero)

    # ---- R-P and R'-P (sector 2@tau + one letter 1, resp. 2@t + one letter 1) ----
    def YP(space, other, flip, sother):
        Dl = lambda o: 2 * L - sother - 2 * s(o)
        gam = osum(space, None, lambda o: s(o) / Dl(o))
        gamt = osum(space, None, lambda o: eps(flip, o) * s(o) / Dl(o))
        rat = gamt / (2 + gam)

        def fA(oc, orr, same):
            if same and oc in space:
                return s(oc) * (kf(Fr(1, 2)) - (1 if oc == flip else 0))
            return zero

        def fB(oc, orr, same):
            if same:
                return s(oc) * rat / 2
            nu = 1 if orr == other else 2
            return s(orr) * s(oc) * nu * (eps(flip, oc) - rat) / (2 * Dl(oc))
        return op(Kk, space, fA), op(Kk, space, fB)

    A_RP, B_RP = YP(Ut, 'T', 't', b)
    A_RpP, B_RpP = YP(UT, 't', 'T', a)

    # ---- Q-Q ----
    fQ = lambda o: LP - 2 * s(o)
    CQ = lambda x, z: LP - 2 * s(x) - 2 * s(z)
    NQ = Nop(U, CQ, None, None)
    Y_QQ = rinv(D2[U] + NQ.scale(Fr(1, 2))).scale(2)
    A_QQ = D2[U].scale(2)
    B_QQ = A_QQ - Y_QQ

    # ---- R-Q and R'-Q ----
    def YQ(space, flip, sother):
        Cp = lambda x, z: 2 * L - sother - 2 * s(x) - 2 * s(z)
        Ah = D2[space] + Nop(space, Cp, None, None).scale(Fr(1, 2))
        H = comp(rinv(Ah), Nop(space, Cp, None, flip)).scale(Fr(1, 2)).restrict(dom=U)
        Aq = op(U, space, lambda oy, ox, same: -s(oy) if same else zero)
        return Aq, H.scale(-1)      # entry = A + (full - A): full = A + H, B-part = -H

    A_RQ, B_RQ = YQ(Ut, 't', b)
    A_RpQ, B_RpQ = YQ(UT, 'T', a)

    # ---- P-Q ----
    gQ = osum(U, None, lambda o: s(o) / fQ(o))
    sQ = 2 + gQ

    def fPQA(orr, oc, same):
        return -s(oc) if same else zero

    def fPQB(orr, oc, same):
        if same:
            return -s(oc) * gQ / sQ
        nu = 1 if orr in ('t', 'T') else 2
        return -2 * nu * s(orr) * s(oc) / (sQ * fQ(oc))
    A_PQ, B_PQ = op(U, Kk, fPQA), op(U, Kk, fPQB)

    def Y(A, B, sc=None):
        sc = bs if sc is None else sc
        return A - B.scale(sc)

    blocks = dict(
        RR=D_R, RpRp=D_Rp, RRp=Y(A_RRp, B_RRp), PP=Y(A_PP, B_PP), RP=Y(A_RP, B_RP), RpP=Y(A_RpP, B_RpP),
        QQ=Y(A_QQ, B_QQ), RQ=Y(A_RQ, B_RQ), RpQ=Y(A_RpQ, B_RpQ), PQ=Y(A_PQ, B_PQ))
    return dict(K=K, Op=Op, op=op, comp=comp, rinv=rinv, blocks=blocks, Ninv=Ninv, spaces=dict(K=Kk, Ut=Ut, UT=UT, U=U),
                G=G, mm=mm, L=L, LV=LV, size=size, sq=sq)


def gram_form(alg, blocks, doms):
    """blocks[i][j]: Op from doms[j] to doms[i].  Returns the invariant form (Gram-weighted reduced matrix) and the
    standard forms (one matrix per level)."""
    mm, G = alg['mm'], alg['G']
    rows = []
    for i, di in enumerate(doms):
        for r in range(len(di)):
            row = []
            for j, dj in enumerate(doms):
                row += mm(G(di), blocks[i][j].M)[r]
            rows.append(row)
    std = [[[blocks[i][j].sd[lv] for j in range(len(doms))] for i in range(len(doms))] for lv in alg['LV']]
    return rows, std


def forms(alg, which='Phi'):
    B = alg['blocks']
    comp = alg['comp']
    sp = alg['spaces']
    if which == 'Phi':
        RP, RpP, PQ, Ni = B['RP'], B['RpP'], B['PQ'], alg['Ninv']
        cols = [RP.T(), RpP.T(), PQ]           # P -> (R, R', Q) couplings, as ops from the three spaces to P
        diag = [[B['RR'], B['RRp'], B['RQ']], [B['RRp'].T(), B['RpRp'], B['RpQ']], [B['RQ'].T(), B['RpQ'].T(), B['QQ']]]
        F = [[diag[i][j] - comp(cols[i].T(), Ni, cols[j]) for j in range(3)] for i in range(3)]
        return gram_form(alg, F, [sp['Ut'], sp['UT'], sp['U']])
    if which == 'G':
        Kn = tuple(o for o in sp['K'] if o != 't')
        RPg = B['RP'].restrict(dom=Kn)
        RpPg = B['RpP'].restrict(dom=Kn)
        PPg = B['PP'].restrict(dom=Kn, cod=Kn)
        PQg = B['PQ'].restrict(cod=Kn)
        F = [[B['RR'], B['RRp'], RPg, B['RQ']], [B['RRp'].T(), B['RpRp'], RpPg, B['RpQ']],
             [RPg.T(), RpPg.T(), PPg, PQg], [B['RQ'].T(), B['RpQ'].T(), PQg.T(), B['QQ']]]
        return gram_form(alg, F, [sp['Ut'], sp['UT'], Kn, sp['U']])
    raise ValueError(which)


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
                A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return d


def symbolic(l, fixed=None, offsets=None, bscale=1, tw=(None, None)):
    """the algebra over Q(a, b, q_1..q_{l-1}, x_i): q_0 = 1, u_i = offsets[i] + x_i, or u_i = fixed[i] exactly.
    tw = (ta, tb): the squared weights of t and tau; None = a free variable ('a', 'b'), or a level number i (the copy
    has the squared weight q_i of level i: two weight values in all)."""
    fixed = fixed or [None] * l
    offsets = offsets or [2] * l
    names = [nm for nm, v in zip(('a', 'b'), tw) if v is None] + [f'q{i}' for i in range(1, l)] +         [f'x{i}' for i in range(l) if fixed[i] is None]
    K = Field(names)
    g = dict(zip(names, K.gens))
    qs = [K(1)] + [g[f'q{i}'] for i in range(1, l)]
    us = [K(fixed[i]) if fixed[i] is not None else g[f'x{i}'] + offsets[i] for i in range(l)]
    a = g['a'] if tw[0] is None else qs[tw[0]]
    b = g['b'] if tw[1] is None else qs[tw[1]]
    return make(K, a, b, qs, us, bscale), K, names


class _FrK:
    """the field Q (Fractions) in place of the symbolic field: the algebra at a rational point."""
    def __call__(self, c):
        if isinstance(c, flint.fmpq):
            return Fr(int(c.p), int(c.q))
        return Fr(c)


def point(a, b, qs, us, bscale=1):
    """the algebra at a rational point (exact, Fractions)."""
    return make(_FrK(), Fr(a), Fr(b), [Fr(x) for x in qs], [Fr(x) for x in us], bscale)


def leading_minors(A, one, callback=None):
    """all leading principal minors D_1, ..., D_n of the symmetric matrix A in one pass of symmetric Gaussian
    elimination without pivoting (D_m = D_{m-1} * pivot_m).  Stops (returns the minors so far) at a zero pivot.
    callback(m, D_m) is called after each minor; if it returns False the elimination stops."""
    A = [list(r) for r in A]
    n = len(A)
    out = []
    D = one
    for m in range(n):
        p = A[m][m]
        if p == 0:
            break
        D = D * p
        out.append(D)
        if callback is not None and callback(m + 1, D) is False:
            break
        if m + 1 == n:
            break
        col = [A[i][m] for i in range(n)]
        for i in range(m + 1, n):
            if col[i] == 0:
                continue
            f = col[i] / p
            for j in range(i, n):
                if col[j] != 0:
                    A[i][j] = A[i][j] - f * col[j]
                    if j != i:
                        A[j][i] = A[i][j]
    return out
