"""Two-level U-weights, symbolic in the multiplicities: the same reductions as two_level.py (invariant part) and
two_level_std.py (standard part of level 1), implemented with SymPy's sparse polynomial rings over the rational
function field Q(w_t, om_1, om_2, u_1, u_2) instead of expression trees (much faster).  w_t, om_i are the square roots
of sigma = s_t and s_i; at the end every coefficient is checked to be even in w_t, om_1, om_2 and converted to a
rational function of (sigma, s1, s2, u2, u1) (python-flint), then the Schur complement is taken by exact elimination.

The closed form used is that of cf.py with x_t = 0 (see equal_u.py for the formula).  Usage:
  python two_level_fast.py inv one     (u1 = 1, u2 >= 2)        -> two_level_one.pkl  (Q3)
  python two_level_fast.py inv gen     (u1, u2 >= 2)            -> two_level_gen.pkl  (Q3)
  python two_level_fast.py std gen     (u1 >= 3, u2 >= 1)       -> two_level_std_gen.pkl  (lam1)
  python two_level_fast.py std u1is2   (u1 = 2, u2 >= 1)        -> two_level_std_u1is2.pkl  (lam1)
Q3 = G_z on the invariant letter-2 data (rho, z1, z2) (sums over the levels); lam1 = G_z[h,h] - G_z[h,h'] for
h != h' in U1.  Invariance/equivariance: the eliminated (C_u) block is positive definite (Theorem C) and invariant
under S_{u1} x S_{u2}, so the minimiser for invariant (resp. level-1 standard) data lies in the invariant (resp.
level-1 standard isotypic) part, which the parametrisations below cover (see two_level.py, two_level_std.py)."""
import sys, os, pickle, time
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sympy.polys.domains import QQ
from sympy.polys.fields import field
from sympy.polys.rings import ring
import flint
from flintfield import Field, F as FF

Kf, WT, O1, O2, U1s, U2s = field('wt o1 o2 u1 u2', QQ)
HALF = Kf.one / 2
OM = {1: O1, 2: O2}
SS = {1: O1 ** 2, 2: O2 ** 2}
SIG = WT ** 2


class Ctx:
    def __init__(self, same, lr, lrp):
        self.same, self.lr, self.lrp = same, lr, (lr if same else lrp)


def make(kind, case):
    if kind == 'inv':
        R, *gens = ring('p1 p2 o11 o12 o21 o22 q rho z1 z2', Kf)
        fr = frp = F2 = None
    else:
        R, *gens = ring('al be ga de ep z fr frp F2', Kf)
        fr, frp, F2 = gens[6], gens[7], gens[8]
    UU = {1: (Kf.one if (kind == 'inv' and case == 'one') else (2 * Kf.one if case == 'u1is2' else U1s)), 2: U2s}
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
        if fr is None:
            return e
        if cx.same:
            e = e.compose(frp, fr)
        if cx.lr == 2:
            e = e.compose(fr, zero)
        if cx.lrp == 2:
            e = e.compose(frp, zero)
        return e

    def merge(v, cx):
        d = {}
        for k, c in v.items():
            kk = 'er' if (cx.same and k == 'erp') else k
            d[kk] = d.get(kk, zero) + fix(R(c), cx)
        return d

    def level_of(k, cx):
        return {'er': cx.lr, 'erp': cx.lrp}.get(k)

    def fval(k, cx):
        return fix({'er': fr, 'erp': frp}[k], cx)

    def wdot(x, y, cx, diagw):
        tot = zero
        for k1, c1 in x.items():
            if k1 == 't':
                continue
            for k2, c2 in y.items():
                if k2 == 't':
                    continue
                a, b = sorted((k1, k2))
                if a.startswith('one') and b.startswith('one'):
                    i, j = int(a[3]), int(b[3])
                    val = UU[i] * diagw(i) if i == j else None
                elif a == 'f' and b == 'f':
                    val = F2 * diagw(1)
                elif (a.startswith('one') and b == 'f') or (b.startswith('one') and a == 'f'):
                    val = None
                elif a.startswith('one') or b.startswith('one'):
                    ko, ke = (a, b) if a.startswith('one') else (b, a)
                    i = int(ko[3])
                    val = diagw(i) if level_of(ke, cx) == i else None
                elif a == 'f' or b == 'f':
                    ke = b if a == 'f' else a
                    val = fval(ke, cx) * diagw(1)
                else:
                    val = diagw(level_of(a, cx)) if a == b else None
                if val is not None:
                    tot += c1 * c2 * val
        return tot

    def lin(x, cx, coef):
        tot = zero
        for k, c in x.items():
            if k == 't' or k == 'f':
                continue
            if k.startswith('one'):
                i = int(k[3])
                tot += c * (UU[i] * coef(i))
            else:
                tot += c * coef(level_of(k, cx))
        return tot

    def coord(x, key, cx):
        x = merge(x, cx)
        key = 'er' if (cx.same and key == 'erp') else key
        l = level_of(key, cx)
        out = x.get(key, zero) + x.get(f'one{l}', zero)
        if 'f' in x:
            out += x['f'] * fval(key, cx)
        return out

    def ipE(x, y, cx, mt, mU, lam, nt, nU):
        x, y = merge(x, cx), merge(y, cx)
        val = x.get('t', zero) * y.get('t', zero) * mt + wdot(x, y, cx, mU)
        nx = x.get('t', zero) * nt + lin(x, cx, nU)
        ny = y.get('t', zero) * nt + lin(y, cx, nU)
        return val - nx * ny * lam

    def red(e, cx):
        """sum over the copies of the context (r = r' in U_i, or r in U_i, r' in U_j, r != r')."""
        e = fix(e, cx)
        i, j = cx.lr, cx.lrp
        if fr is None:
            if cx.same:
                return e * UU[i]
            cnt = UU[i] * (UU[j] - (1 if i == j else 0))
            return e * cnt
        nf = R.ngens
        out = zero
        for mon, co in e.terms():
            a, b, c = mon[6], mon[7], mon[8]
            rest = R({tuple(mon[:6]) + (0, 0, 0): co})
            deg = a + b + 2 * c
            if deg == 1:
                continue
            if deg not in (0, 2):
                raise ValueError('degree in f > 2')
            if cx.same:
                if c == 1:
                    val = F2 * UU[i]
                elif (a, b) == (0, 0):
                    val = R(UU[i])
                elif (a, b) == (2, 0):
                    val = F2 if i == 1 else None
                else:
                    raise ValueError('unexpected monomial in the same-copy case')
            else:
                cnt = UU[i] * (UU[j] - (1 if i == j else 0))
                if c == 1:
                    val = F2 * cnt
                elif (a, b) == (0, 0):
                    val = R(cnt)
                elif (a, b) == (2, 0):
                    val = F2 * (UU[j] - (1 if i == j else 0)) if i == 1 else None
                elif (a, b) == (0, 2):
                    val = F2 * (UU[i] - (1 if i == j else 0)) if j == 1 else None
                elif (a, b) == (1, 1):
                    val = -F2 if (i == 1 and j == 1) else None
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
    K = SIG + S / 2
    d = {i: SS[i] + SIG + S / 2 for i in (1, 2)}
    lamT = Kf.one / (2 + sum((UU[i] * SS[i] / d[i] for i in (1, 2)), Kf.zero))
    d0t = SIG + (SIG + S) / 2
    d0 = {i: SS[i] + (SIG + S) / 2 for i in (1, 2)}
    lam0 = Kf.one / (2 + SIG / d0t + sum((UU[i] * SS[i] / d0[i] for i in (1, 2)), Kf.zero))

    def Tp(psi, key, cx):
        l = level_of('er' if (cx.same and key == 'erp') else key, cx)
        wpsi = lin(merge(psi, cx), cx, lambda i: OM[i])
        c = coord(psi, key, cx)
        e = {key: R.one}
        return add(mul(wpsi * HALF, e), mul(R(OM[l] / 2), psi), mul(c * OM[l], e))

    if kind == 'inv':
        p1, p2, o11, o12, o21, o22, q, rho, z1, z2 = gens
        P = {1: p1, 2: p2}
        O = {(1, 1): (zero if UU[1] == Kf.one else o11), (1, 2): o12, (2, 1): o21, (2, 2): o22}

        def Xs(i, j, same):
            return (P[i] if (same and i == j) else O[(i, j)]) * (OM[i] * OM[j])

        def row(i, key):
            v = add({'one1': Xs(i, 1, False)}, {'one2': Xs(i, 2, False)})
            return add(v, {key: Xs(i, i, True) - Xs(i, i, False)})

        A = sum_rr(lambda cx: Xs(cx.lr, cx.lrp, cx.same) ** 2) * HALF + sum_r(lambda cx: Xs(cx.lr, cx.lr, True) ** 2)
        trX = sum((Xs(i, i, True) * UU[i] for i in (1, 2)), zero)
        wXw = zero
        for i in (1, 2):
            wXw += Xs(i, i, True) * (UU[i] * SS[i])
            for j in (1, 2):
                cnt = UU[i] * (UU[j] - (1 if i == j else 0))
                wXw += Xs(i, j, False) * (cnt * OM[i] * OM[j])
        zeta = {'one1': z1 * (WT * O1), 'one2': z2 * (WT * O2)}
        wzeta = z1 * (UU[1] * O1 * WT * O1) + z2 * (UU[2] * O2 * WT * O2)
        zz = z1 ** 2 * (UU[1] * SIG * SS[1]) + z2 ** 2 * (UU[2] * SIG * SS[2])
        A += (trX + rho) ** 2 * HALF + rho ** 2 * (3 * HALF) + zz
        cpl = q * 2 * (-trX + wXw * (Kf.one / (2 * K)) - (rho * (2 * K) + wzeta * WT) * (Kf.one / (2 * K - SIG))) + q ** 2 * 2
        y_extra = {'t': rho}
        wy = rho * WT + wzeta
        elim = [p1, p2, o12, o21, o22, q] + ([] if UU[1] == Kf.one else [o11])
        keep = [rho, z1, z2]
    else:
        al, be, ga, de, ep, z = gens[:6]
        if case == 'u1is2':
            ga = zero

        def Xs(i, j, same, a, b):
            pre = OM[i] * OM[j]
            if i == 1 and j == 1:
                return ((al + be + ga) * a if same else be * a + ga * b) * pre
            if i == 1 and j == 2:
                return de * a * pre
            if i == 2 and j == 1:
                return ep * b * pre
            return zero

        def row(i, key):
            a = fr if key == 'er' else frp
            if i == 1:
                return {'one1': be * a * SS[1], 'f': ga * SS[1], 'one2': de * a * (O1 * O2), key: al * a * SS[1]}
            return {'f': ep * (O2 * O1)}

        A = sum_rr(lambda cx: Xs(cx.lr, cx.lrp, cx.same, fr, (fr if cx.same else frp)) ** 2) * HALF
        A += sum_r(lambda cx: Xs(cx.lr, cx.lr, True, fr, fr) ** 2)
        zeta = {'f': z * (WT * O1)}
        A += z ** 2 * F2 * (SIG * SS[1])
        cpl = zero
        y_extra = {}
        wy = zero
        elim = [al, be, de, ep] if case == 'u1is2' else [al, be, ga, de, ep]
        keep = [z]

    def BU(cx):
        psi_r = row(cx.lr, 'er')
        psi_rp = row(cx.lrp, 'erp') if not cx.same else row(cx.lr, 'er')
        a1 = Tp(psi_r, 'erp', cx)
        a2 = Tp(psi_rp, 'er', cx)
        return ipE(a1, a2, cx, Kf.zero, lambda i: Kf.one / d[i], lamT, Kf.zero, lambda i: OM[i] / d[i])

    def cross(cx):
        i = cx.lr
        psi_r = row(i, 'er')
        wpsi = lin(merge(psi_r, cx), cx, lambda l: OM[l])
        T2 = add(mul(wpsi * HALF, {'t': R.one}), mul(R(-WT / 2), psi_r))
        y = add(y_extra, zeta)
        yr = coord(zeta, 'er', cx)
        er = {'er': R.one}
        T1y = add(mul(wy * HALF, er), mul(R(OM[i] / 2), y), mul(yr * OM[i], er))
        return ipE(T2, T1y, cx, Kf.one / d0t, lambda l: Kf.one / d0[l], lam0, WT / d0t, lambda l: OM[l] / d0[l])

    Phi = A + cpl - sum_rr(BU) - sum_r(cross) * 2
    return R, Phi, elim, keep, F2


def coeff_matrix(R, Phi, allv, F2):
    """Hessian/2 of the quadratic form Phi in the variables allv (divided by F2 for the standard part)."""
    idx = {R.gens.index(v): k for k, v in enumerate(allv)}
    n = len(allv)
    H = [[Kf.zero] * n for _ in range(n)]
    fpos = R.gens.index(F2) if F2 is not None else None
    for mon, co in Phi.terms():
        if fpos is not None:
            if mon[fpos] != 1:
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


def to_flint(c, K, case_vars):
    """c in Kf (variables wt, o1, o2, u1, u2): halve the exponents of wt, o1, o2 (they must be even)."""
    def conv(p):
        out = {}
        for mon, co in p.terms():
            if mon[0] % 2 or mon[1] % 2 or mon[2] % 2:
                raise ValueError('odd exponent of a square root')
            m = {'sigma': mon[0] // 2, 's1': mon[1] // 2, 's2': mon[2] // 2, 'u1': mon[3], 'u2': mon[4]}
            key = tuple(m[v] for v in case_vars)
            if any(m[v] for v in m if v not in case_vars):
                raise ValueError('variable not expected')
            fq = Fraction(str(co))
            out[key] = flint.fmpq(fq.numerator, fq.denominator)
        return K.ctx.from_dict(out)
    return FF(K, conv(c.numer), conv(c.denom))


def run(kind, case):
    t0 = time.time()
    R, Phi, elim, keep, F2 = make(kind, case)
    print(f'{kind} {case}: form built ({time.time() - t0:.1f} s, {len(Phi.terms())} terms)', flush=True)
    allv = elim + keep
    Hk = coeff_matrix(R, Phi, allv, F2)
    case_vars = ['sigma', 's1', 's2', 'u2'] + ([] if (case in ('one', 'u1is2')) else ['u1'])
    K = Field(case_vars)
    n = len(allv)
    H = [[to_flint(Hk[i][j], K, case_vars) for j in range(n)] for i in range(n)]
    print(f'  converted ({time.time() - t0:.1f} s); eliminating {len(elim)}', flush=True)
    ne = len(elim)
    for p in range(ne):
        pv = H[p][p]
        assert pv != 0
        for i in range(p + 1, n):
            if H[i][p] == 0:
                continue
            f = H[i][p] / pv
            for j in range(p + 1, n):
                if H[p][j] != 0:
                    H[i][j] = H[i][j] - f * H[p][j]
    nk = len(keep)
    res = [[H[ne + i][ne + j] for j in range(nk)] for i in range(nk)]
    sym = all(res[i][j] == res[j][i] for i in range(nk) for j in range(nk))
    out = [[(res[i][j].n.to_dict(), res[i][j].d.to_dict()) for j in range(nk)] for i in range(nk)]
    if kind == 'inv':
        fn = f'two_level_{case}.pkl'
        pickle.dump(dict(Q3=out, names=case_vars, case=case), open(fn, 'wb'))
    else:
        fn = f'two_level_std_{case}.pkl'
        pickle.dump(dict(lam1=out[0][0], names=case_vars, case=case), open(fn, 'wb'))
    print(f'  done ({time.time() - t0:.1f} s); symmetric {sym}; wrote {fn}', flush=True)
    return sym


if __name__ == '__main__':
    ok = run(sys.argv[1], sys.argv[2])
    print('PASS' if ok else 'FAIL')
