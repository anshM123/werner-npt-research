"""Two U-weight levels and two copies carrying the letter 2 (plus further T-copies with zero data): the blocks G_tt and
G_tt' of the Schur complement G, symbolically in the multiplicities.

Setting (relative-anchors/README.md Section 3; typeI_G_em/NOTES.md Section 1).  U = U1 + U2, |U_i| = u_i, every copy of
U_i has squared weight s_i = om_i^2; by homogeneity (2^Z Y does not change when all weights are multiplied by one
positive factor) s2 = 1.  T contains the copies A, B (squared weights a = wa^2, b = wb^2) and further copies with zero
letter-2 data; sigma = s_T is a free symbol (sigma >= a + b).  S = u1 s1 + u2 s2, K = sigma + S/2.
In the gauge sum_t w_t psi_t = 0 (typeI_G_em/NOTES.md Section 1):
    G_tt  = diag(2, I) - (1/4)(1 - s_t/sigma) L_t^T L_t - C_t^T C^{-1} C_t,
    G_tt' = (1/2) e_rho e_rho^T - P_tt' + (w_t w_t'/(4 sigma)) L_t^T L_t' - C_t^T C^{-1} C_t'      (t != t').
The quadratic form below is the sum of all these terms before the elimination of (X_U, q):
    Phi = (C_u form on (X_U, q)) + sum_t (coupling of y_t to (X_U, q)) + (letter-2 terms),
and the Schur complement onto the letter-2 coordinates gives the 2 x 2 array of blocks (t, t') for t, t' in {A, B}.

The block is invariant under S_{u1} x S_{u2} (permutations of each level), C is positive definite (Theorem C), so the
minimiser over (X_U, q) is unique and equivariant.  The letter-2 data of each copy t split into
  * the invariant part: rho_t, zeta_t = w_t (om_1 z_t^1 1_{U1} + om_2 z_t^2 1_{U2}); minimiser x_hh = p_i (h in U_i),
    x_rh = o_ij (r in U_i, h in U_j, r != h; o_11 absent if u1 = 1), q            (kind 'inv');
  * the standard part of level 1: rho = 0, zeta_t = w_t om_1 z_t f, f on U1 with sum f = 0; minimiser
    x = al diag f + be f 1^T + ga 1 f^T on U1 x U1 (ga absent if u1 = 2), x_rh = de f_r on U1 x U2, x_hr = ep f_r on
    U2 x U1, zero on U2 x U2, q = 0 (the standard irrep of S_{u1} occurs 3 (resp. 2) + 1 + 1 times in R^{U x U} and not
    in q); value per |f|^2                                                          (kind 'std1').
  (The standard part of level 2 is the same computation with the levels exchanged; by homogeneity it is covered by
  'std1' with s1 -> 1/s1, see NOTES.md.)
Scaled coordinates X_rh = w_r w_h x_rh, zeta_{t,h} = w_t w_h z_{t,h}.

Usage (append eqT for all T-weights equal to s2 = 1, i.e. a = b = 1 and sigma = m):
        python tl_m.py inv one    (u1 = 1, u2 >= 2 symbolic)       -> tl_inv_one.pkl
        python tl_m.py inv gen    (u1, u2 >= 2 symbolic)           -> tl_inv_gen.pkl
        python tl_m.py std1 gen   (u1 >= 3 symbolic, u2 >= 1)      -> tl_std1_gen.pkl
        python tl_m.py std1 two   (u1 = 2, u2 >= 1 symbolic)       -> tl_std1_two.pkl
Output: G as a matrix of flint rational functions (serialised) in the variables NAMES = (u1, u2, sig, s1, a, b) (u1
absent when fixed), order inv: (rhoA, zA1, zA2, rhoB, zB1, zB2); std1: (zA, zB)."""
import sys, os, pickle, time
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sympy.polys.domains import QQ
from sympy.polys.fields import field
from sympy.polys.rings import ring
import flint
from flintfield import Field, F as FF

Kf, WA, WB, O1, U1s, U2s, SIG = field('wa wb o1 u1 u2 sig', QQ)
O2 = Kf.one                      # s2 = 1
HALF = Kf.one / 2
OM = {1: O1, 2: O2}
SS = {1: O1 ** 2, 2: O2 ** 2}
# T-weight modes: general (a, b symbolic); eqT = T22 (a = b = 1 = s2); T11 (a = b = s1); T12 (a = s1, b = 1)
TMODE = next((x for x in sys.argv[3:] if x in ('eqT', 'T11', 'T12', 'T22')), None)
EQT = TMODE is not None
WT = {None: {'A': WA, 'B': WB}, 'eqT': {'A': Kf.one, 'B': Kf.one}, 'T22': {'A': Kf.one, 'B': Kf.one},
      'T11': {'A': O1, 'B': O1}, 'T12': {'A': O1, 'B': Kf.one}}[TMODE]
STW = {t: WT[t] ** 2 for t in WT}
TK = ('A', 'B')


# (u1, u2): an integer, or None for symbolic
CASES = {'inv': dict([('one', (1, None)), ('gen', (None, None))] + [('f%d' % n, (n, None)) for n in range(2, 61)]),
         'std1': {'gen': (None, None), 'two': (2, None)},
         'std2': {'one': (1, None), 'onetwo': (1, 2), 'gen': (None, None), 'two': (None, 2)}}


class Ctx:
    def __init__(self, same, lr, lrp):
        self.same, self.lr, self.lrp = same, lr, (lr if same else lrp)


def make(kind, case):
    LV = {'inv': None, 'std1': 1, 'std2': 2}[kind]      # level carrying the standard vector f
    OV = None if LV is None else 3 - LV
    u1c, u2c = CASES[kind][case]
    UU = {1: (U1s if u1c is None else u1c * Kf.one), 2: (U2s if u2c is None else u2c * Kf.one)}
    if kind == 'inv':
        R, *gens = ring('p1 p2 o11 o12 o21 o22 q rA zA1 zA2 rB zB1 zB2', Kf)
        fr = frp = F2 = None
    else:
        R, *gens = ring('al be ga de ep zA zB fr frp F2', Kf)
        fr, frp, F2 = gens[-3:]
    zero = R.zero

    def add(*vs):
        d = {}
        for v in vs:
            for k, c in v.items():
                d[k] = d.get(k, zero) + c
        return d

    def mul(c, v):
        return {k: c * x for k, x in v.items()}

    def fix(e, cx):
        e = R(e)
        if fr is None:
            return e
        if cx.same:
            e = e.compose(frp, fr)
        if cx.lr != LV:
            e = e.compose(fr, zero)
        if cx.lrp != LV:
            e = e.compose(frp, zero)
        return e

    def merge(v, cx):
        d = {}
        for k, c in v.items():
            kk = 'er' if (cx.same and k == 'erp') else k
            d[kk] = d.get(kk, zero) + fix(c, cx)
        return d

    def lev(k, cx):
        return {'er': cx.lr, 'erp': cx.lrp}[k]

    def fval(k, cx):
        return fix({'er': fr, 'erp': frp}[k], cx)

    def dot(x, y, cx, dU, dT):
        """<x, D y>, D = diag(dU(level) on U, dT(key) on the T-copies); x, y merged."""
        tot = zero
        for k1, c1 in x.items():
            for k2, c2 in y.items():
                if k1 in TK or k2 in TK:
                    if k1 == k2:
                        tot += c1 * c2 * dT(k1)
                    continue
                a, b = sorted((k1, k2))
                if a.startswith('one') and b.startswith('one'):
                    i, j = int(a[3]), int(b[3])
                    val = UU[i] * dU(i) if i == j else None
                elif a == 'f' and b == 'f':
                    val = F2 * dU(LV)
                elif (a.startswith('one') and b == 'f') or (b.startswith('one') and a == 'f'):
                    val = None
                elif a.startswith('one') or b.startswith('one'):
                    ko, ke = (a, b) if a.startswith('one') else (b, a)
                    i = int(ko[3])
                    val = dU(i) if lev(ke, cx) == i else None
                elif a == 'f' or b == 'f':
                    ke = b if a == 'f' else a
                    val = fval(ke, cx) * dU(LV)
                else:
                    val = dU(lev(a, cx)) if a == b else None
                if val is not None:
                    tot += c1 * c2 * val
        return tot

    def lin(x, cx, cU, cT):
        """sum_h c(h) x_h with c = cU(level) on U and cT(key) on the T-copies; x merged."""
        tot = zero
        for k, c in x.items():
            if k == 'f':
                continue
            if k in TK:
                tot += c * cT(k)
            elif k.startswith('one'):
                i = int(k[3])
                tot += c * (UU[i] * cU(i))
            else:
                tot += c * cU(lev(k, cx))
        return tot

    def coord(x, key, cx):
        x = merge(x, cx)
        key = 'er' if (cx.same and key == 'erp') else key
        l = lev(key, cx)
        out = x.get(key, zero) + x.get(f'one{l}', zero)
        if 'f' in x:
            out += x['f'] * fval(key, cx)
        return out

    class Emat:
        """E = diag(d) + 1/2 wt wt^T on U + (T-copies in `tk`), d_x = s_x + K - shift."""
        def __init__(self, tk, shift):
            self.tk = tk
            self.dU = {i: SS[i] + Kc - shift for i in (1, 2)}
            self.dT = {t: STW[t] + Kc - shift for t in tk}
            g = sum((UU[i] * SS[i] / self.dU[i] for i in (1, 2)), Kf.zero)
            g += sum((STW[t] / self.dT[t] for t in tk), Kf.zero)
            self.lam = Kf.one / (2 + g)

        def ip(self, x, y, cx):
            x, y = merge(x, cx), merge(y, cx)
            for k in list(x) + list(y):
                if k in TK and k not in self.tk:
                    raise ValueError('vector component outside the free set')
            val = dot(x, y, cx, lambda i: Kf.one / self.dU[i], lambda t: Kf.one / self.dT[t])
            ex = lin(x, cx, lambda i: OM[i] / self.dU[i], lambda t: WT[t] / self.dT[t])
            ey = lin(y, cx, lambda i: OM[i] / self.dU[i], lambda t: WT[t] / self.dT[t])
            return val - ex * ey * self.lam

    def red(e, cx):
        """sum over the copies of the context (r = r' in U_i, or r in U_i, r' in U_j, r != r')."""
        e = fix(e, cx)
        i, j = cx.lr, cx.lrp
        if fr is None:
            if cx.same:
                return e * UU[i]
            return e * (UU[i] * (UU[j] - (1 if i == j else 0)))
        ia, ib, ic = R.gens.index(fr), R.gens.index(frp), R.gens.index(F2)
        out = zero
        for mon, co in e.terms():
            a, b, c = mon[ia], mon[ib], mon[ic]
            mm = list(mon)
            mm[ia] = mm[ib] = mm[ic] = 0
            rest = R({tuple(mm): co})
            deg = a + b + 2 * c
            if deg % 2 == 1:
                if deg == 1:
                    continue
                raise ValueError('odd degree > 1 in f')
            if deg not in (0, 2):
                raise ValueError('degree in f > 2')
            if cx.same:
                if c == 1:
                    val = F2 * UU[i]
                elif (a, b) == (0, 0):
                    val = R(UU[i])
                elif (a, b) == (2, 0):
                    val = F2 if i == LV else None
                else:
                    raise ValueError('unexpected monomial (same)')
            else:
                cnt = UU[i] * (UU[j] - (1 if i == j else 0))
                if c == 1:
                    val = F2 * cnt
                elif (a, b) == (0, 0):
                    val = R(cnt)
                elif (a, b) == (2, 0):
                    val = F2 * (UU[j] - (1 if i == j else 0)) if i == LV else None
                elif (a, b) == (0, 2):
                    val = F2 * (UU[i] - (1 if i == j else 0)) if j == LV else None
                elif (a, b) == (1, 1):
                    val = -F2 if (i == LV and j == LV) else None
                else:
                    val = None
            if val is not None:
                out += rest * val
        return out

    def sum_rr(f):
        tot = zero
        for i in (1, 2):
            cx = Ctx(True, i, i)
            tot += red(f(cx), cx)
            for j in (1, 2):
                if i == j and UU[i] == Kf.one:
                    continue
                cx = Ctx(False, i, j)
                tot += red(f(cx), cx)
        return tot

    def sum_r(f):
        tot = zero
        for i in (1, 2):
            cx = Ctx(True, i, i)
            tot += red(f(cx), cx)
        return tot

    S = UU[1] * SS[1] + UU[2] * SS[2]
    Kc = SIG + S / 2
    ET = Emat((), Kf.zero)                                   # pattern T: free set U, d_h = s_h + K
    Et = {t: Emat((t,), STW[t] / 2) for t in TK}              # pattern T - t: free U + t
    Ett = Emat(TK, (STW['A'] + STW['B']) / 2)                 # pattern T - {A, B}

    def Tp(psi, key, cx):
        """T_{(r',1)} on a pattern-T vector psi (on U), r' the copy of `key`."""
        l = lev('er' if (cx.same and key == 'erp') else key, cx)
        wpsi = lin(merge(psi, cx), cx, lambda i: OM[i], None)
        c = coord(psi, key, cx)
        e = {key: R.one}
        return add(mul(wpsi * HALF, e), mul(R(OM[l] / 2), psi), mul(c * OM[l], e))

    # ---------------- parametrisation
    if kind == 'inv':
        p1, p2, o11, o12, o21, o22, q, rA, zA1, zA2, rB, zB1, zB2 = gens
        P = {1: p1, 2: p2}
        O = {(1, 1): (zero if UU[1] == Kf.one else o11), (1, 2): o12, (2, 1): o21, (2, 2): o22}

        def Xs(i, j, same):
            return (P[i] if (same and i == j) else O[(i, j)]) * (OM[i] * OM[j])

        def row(i, key):
            v = add({'one1': Xs(i, 1, False)}, {'one2': Xs(i, 2, False)})
            return add(v, {key: Xs(i, i, True) - Xs(i, i, False)})

        Apart = sum_rr(lambda cx: Xs(cx.lr, cx.lrp, cx.same) ** 2) * HALF + sum_r(lambda cx: Xs(cx.lr, cx.lr, True) ** 2)
        trX = sum((Xs(i, i, True) * UU[i] for i in (1, 2)), zero)
        wXw = zero
        for i in (1, 2):
            wXw += Xs(i, i, True) * (UU[i] * SS[i])
            for j in (1, 2):
                cnt = UU[i] * (UU[j] - (1 if i == j else 0))
                wXw += Xs(i, j, False) * (cnt * OM[i] * OM[j])
        Apart += trX ** 2 * HALF
        Cq = q * 2 * (-trX + wXw * (Kf.one / (2 * Kc))) + q ** 2 * 2
        rho = {'A': rA, 'B': rB}
        cz = {'A': {1: zA1 * (WT['A'] * O1), 2: zA2 * (WT['A'] * O2)}, 'B': {1: zB1 * (WT['B'] * O1), 2: zB2 * (WT['B'] * O2)}}
        zeta = {t: {'one1': cz[t][1], 'one2': cz[t][2]} for t in TK}
        elim = [p1, p2, o12, o21, o22, q] + ([] if UU[1] == Kf.one else [o11])
        keep = [rA, zA1, zA2, rB, zB1, zB2]
    else:
        al, be, ga, de, ep, zA, zB = gens[:7]
        gamma_absent = UU[LV] == 2 * Kf.one
        if gamma_absent:
            ga = zero

        def Xs(i, j, same, fa, fb):
            pre = OM[i] * OM[j]
            if i == LV and j == LV:
                return ((al + be + ga) * fa if same else be * fa + ga * fb) * pre
            if i == LV and j == OV:
                return de * fa * pre
            if i == OV and j == LV:
                return ep * fb * pre
            return zero

        def row(i, key):
            fa = fr if key == 'er' else frp
            if i == LV:
                return {f'one{LV}': be * fa * SS[LV], 'f': ga * SS[LV], f'one{OV}': de * fa * (OM[LV] * OM[OV]),
                        key: al * fa * SS[LV]}
            return {'f': ep * (OM[OV] * OM[LV])}

        Apart = sum_rr(lambda cx: Xs(cx.lr, cx.lrp, cx.same, fr, (fr if cx.same else frp)) ** 2) * HALF
        Apart += sum_r(lambda cx: Xs(cx.lr, cx.lr, True, fr, fr) ** 2)
        trX = zero
        Cq = zero
        rho = {'A': zero, 'B': zero}
        cz = {'A': {LV: zA * (WT['A'] * OM[LV])}, 'B': {LV: zB * (WT['B'] * OM[LV])}}
        zeta = {t: {'f': cz[t][LV]} for t in TK}
        elim = [al, be, de, ep] if gamma_absent else [al, be, ga, de, ep]
        keep = [zA, zB]

    # ---------------- the (C_u) B-term
    def BU(cx):
        psi_r = row(cx.lr, 'er')
        psi_rp = row(cx.lrp, 'erp') if not cx.same else row(cx.lr, 'er')
        return ET.ip(Tp(psi_r, 'erp', cx), Tp(psi_rp, 'er', cx), cx)

    PhiC = Apart + Cq - sum_rr(BU)

    # ---------------- couplings of y_t to (X_U, q)
    def wU_dot(v, cx):
        return lin(merge(v, cx), cx, lambda i: OM[i], None)

    def coupling(t):
        wt, st = WT[t], STW[t]
        cxi = Ctx(True, 1, 1)
        wz = lin(zeta[t], cxi, lambda i: OM[i], None) if kind == 'inv' else zero
        out = trX * rho[t]
        if kind == 'inv':
            out += -2 * q * (2 * Kc * rho[t] + wt * wz) * (Kf.one / (2 * Kc - st))

        def cross(cx):
            i = cx.lr
            psi_r = row(i, 'er')
            T2 = add(mul(wU_dot(psi_r, cx) * HALF, {t: R.one}), mul(R(-wt / 2), psi_r))
            y = add({t: rho[t]}, zeta[t])
            wy = rho[t] * wt + wU_dot(zeta[t], cx)
            yr = coord(zeta[t], 'er', cx)
            er = {'er': R.one}
            T1y = add(mul(wy * HALF, er), mul(R(OM[i] / 2), y), mul(yr * OM[i], er))
            return Et[t].ip(T2, T1y, cx)
        return out - 2 * sum_r(cross)

    # ---------------- letter-2 terms: direct, psi_t elimination (ell), pair
    def ell(t):
        """the coefficients of l_t on the level vectors (inv: {1, 2}) or on f (std: {1})."""
        wt, st = WT[t], STW[t]
        E = Et[t]
        if kind == 'inv':
            c = cz[t]
            vt = wt * rho[t] + sum((UU[i] * OM[i] * c[i] for i in (1, 2)), zero) * HALF
            vU = {i: -wt * c[i] * HALF for i in (1, 2)}
            ev = vt * (wt / E.dT[t]) + sum((UU[i] * (OM[i] / E.dU[i]) * vU[i] for i in (1, 2)), zero)
            Etv = vt * (Kf.one / E.dT[t]) - ev * (wt / E.dT[t] * E.lam)
            EUv = {i: vU[i] * (Kf.one / E.dU[i]) - ev * (OM[i] / E.dU[i] * E.lam) for i in (1, 2)}
            return {i: -c[i] - Etv * OM[i] + wt * EUv[i] for i in (1, 2)}
        return {LV: -cz[t][LV] * (1 + st / (2 * E.dU[LV]))}

    def norm2(lv):
        if kind == 'inv':
            return sum((UU[i] * lv[i] ** 2 for i in (1, 2)), zero)
        return lv[LV] ** 2 * F2

    L = {t: ell(t) for t in TK}
    if kind == 'inv':
        direct = 2 * rA ** 2 + 2 * rB ** 2 + rA * rB
        direct += sum((norm2(cz[t]) for t in TK), zero)
    else:
        direct = sum((cz[t][LV] ** 2 for t in TK), zero) * F2
    psi_part = -sum((norm2(L[t]) for t in TK), zero) / 4
    comb = {i: sum((WT[t] * L[t][i] for t in TK), zero) for i in L['A']}
    psi_part += norm2(comb) * (Kf.one / (4 * SIG))

    def pair():
        """<T_(B,2) y_A, E_{T-AB}^{-1} T_(A,2) y_B>"""
        cx = Ctx(True, 1, 1)
        yA = add({'A': rho['A']}, zeta['A'])
        yB = add({'B': rho['B']}, zeta['B'])
        wA = rho['A'] * WT['A'] + wU_dot(zeta['A'], cx)
        wB = rho['B'] * WT['B'] + wU_dot(zeta['B'], cx)
        v1 = add(mul(wA * HALF, {'B': R.one}), mul(R(-WT['B'] / 2), yA))
        v2 = add(mul(wB * HALF, {'A': R.one}), mul(R(-WT['A'] / 2), yB))
        val = Ett.ip(v1, v2, cx)
        if kind != 'inv':
            # f-components only (sum f = 0 kills the J-part); dot gives F2 already
            pass
        return val

    Phi = PhiC + coupling('A') + coupling('B') + direct + psi_part - 2 * pair()
    return R, Phi, elim, keep, F2


def coeff_matrix(R, Phi, allv, F2):
    idx = {R.gens.index(v): k for k, v in enumerate(allv)}
    n = len(allv)
    H = [[Kf.zero] * n for _ in range(n)]
    fpos = R.gens.index(F2) if F2 is not None else None
    for mon, co in Phi.terms():
        if fpos is not None and mon[fpos] != 1:
            raise ValueError(f'term not proportional to F2: {mon}')
        vs = []
        for pos, e in enumerate(mon):
            if pos == fpos:
                continue
            if e:
                if pos not in idx:
                    raise ValueError(f'unexpected variable in {mon}')
                vs += [idx[pos]] * e
        if len(vs) != 2:
            raise ValueError(f'not a quadratic form: {mon}')
        a, b = vs
        if a == b:
            H[a][a] += co
        else:
            H[a][b] += co / 2
            H[b][a] += co / 2
    return H


def to_flint(c, K, names):
    """c in Kf (wa, wb, o1, u1, u2, sig): halve the exponents of wa, wb, o1 (they must be even)."""
    def conv(p):
        out = {}
        for mon, co in p.terms():
            if mon[0] % 2 or mon[1] % 2 or mon[2] % 2:
                raise ValueError('odd exponent of a square root')
            m = {'a': mon[0] // 2, 'b': mon[1] // 2, 's1': mon[2] // 2, 'u1': mon[3], 'u2': mon[4], 'sig': mon[5]}
            if any(m[v] for v in m if v not in names):
                raise ValueError('variable not expected: %s' % (m,))
            key = tuple(m[v] for v in names)
            fq = Fraction(str(co))
            out[key] = out.get(key, 0) + flint.fmpq(fq.numerator, fq.denominator)
        return K.ctx.from_dict(out)
    return FF(K, conv(c.numer), conv(c.denom))


def ser_poly(p):
    return {tuple(int(e) for e in m): (int(c.p), int(c.q)) for m, c in zip(p.monoms(), p.coeffs())}


def run(kind, case):
    t0 = time.time()
    R, Phi, elim, keep, F2 = make(kind, case)
    print(f'{kind} {case}: form built ({time.time() - t0:.1f} s, {len(Phi.terms())} terms)', flush=True)
    allv = elim + keep
    Hk = coeff_matrix(R, Phi, allv, F2)
    u1c, u2c = CASES[kind][case]
    names = (['u1'] if u1c is None else []) + (['u2'] if u2c is None else []) + ['sig', 's1'] + ([] if EQT else ['a', 'b'])
    K = Field(names)
    n = len(allv)
    ne = len(elim)
    nk = len(keep)
    half = nk // 2
    # needed entries (upper triangle): everything touching an eliminated variable, G_AA (upper) and G_AB
    need = {(i, j) for i in range(n) for j in range(i, n)
            if i < ne or (i - ne < half)}
    H = {}
    for (i, j) in need:
        H[(i, j)] = to_flint(Hk[i][j], K, names)
    print(f'  converted ({time.time() - t0:.1f} s); eliminating {ne}', flush=True)
    for p in range(ne):
        pv = H[(p, p)]
        assert pv != 0
        for i in range(p + 1, n):
            hpi = H[(p, i)]
            if hpi == 0:
                continue
            f = hpi / pv
            for j in range(i, n):
                if (i, j) not in need:
                    continue
                hpj = H[(p, j)]
                if hpj != 0:
                    H[(i, j)] = H[(i, j)] - f * hpj
        for j in range(p, n):
            H.pop((p, j), None)
        print(f'  pivot {p} done ({time.time() - t0:.1f} s)', flush=True)
    out = {}
    for i in range(half):
        for j in range(i, nk):
            g = H[(ne + i, ne + j)]
            out[(i, j)] = (ser_poly(g.n), ser_poly(g.d))
            print(f'  G[{keep[i]},{keep[j]}]: num {len(g.n)} terms, den {len(g.d)} terms', flush=True)
    fn = f'tl_{kind}_{case}' + (f'_{TMODE}' if EQT else '') + '.pkl'
    pickle.dump(dict(G=out, names=names, kind=kind, case=case, keep=[str(v) for v in keep],
                     elim=[str(v) for v in elim]), open(fn, 'wb'))
    print(f'  done ({time.time() - t0:.1f} s); wrote {fn}', flush=True)
    return True


if __name__ == '__main__':
    ok = run(sys.argv[1], sys.argv[2])
    print('PASS' if ok else 'FAIL')
