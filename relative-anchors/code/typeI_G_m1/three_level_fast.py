"""Three levels: two single U-copies with arbitrary weights s1, s2 and u3 >= 2 further U-copies of a common weight s3
(u3 symbolic).  Invariant part only (S_{u3}-invariant letter-2 data (rho, z1, z2, z3)); same construction as
two_level_fast.py (polynomial rings over Q(w_t, om_1, om_2, om_3, u3), then exact elimination with python-flint).
sigma = s_t is put to 1 (homogeneity, typeI_G_u3/NOTES.md Section 1).  Unknowns of the invariant minimiser: p_i = x_hh (h in level i), o_ij = x_rh (r in level i, h in level j, r != h; o_11
and o_22 do not exist), and q.  Output three_level.pkl with Q4 = G_z summed over the levels:
  Q4[0,0] = G_z[rho,rho], Q4[0,i] = sum_{h in level i} G_z[rho,h], Q4[i,j] = sum_{h in level i, h' in level j} G_z[h,h'].
Usage: python three_level_fast.py"""
import sys, os, pickle, time
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sympy.polys.domains import QQ
from sympy.polys.fields import field
from sympy.polys.rings import ring
import flint
from flintfield import Field, F as FF

Kf, O1, O2, O3, U3s = field('o1 o2 o3 u3', QQ)
WT = Kf.one            # sigma = s_t = 1 (homogeneity: G_z(l s) = D_l G_z(s) D_l, D_l = diag(1, l, ..., l))
HALF = Kf.one / 2
LV = (1, 2, 3)
OM = {1: O1, 2: O2, 3: O3}
SS = {i: OM[i] ** 2 for i in LV}
UU = {1: Kf.one, 2: Kf.one, 3: U3s}
SIG = WT ** 2


class Ctx:
    def __init__(self, same, lr, lrp):
        self.same, self.lr, self.lrp = same, lr, (lr if same else lrp)


def make():
    names = 'p1 p2 p3 o12 o13 o21 o23 o31 o32 o33 q rho z1 z2 z3'
    R, *gens = ring(names, Kf)
    G = dict(zip(names.split(), gens))
    zero = R.zero

    def add(*vs):
        d = {}
        for v in vs:
            for k, c in v.items():
                d[k] = d.get(k, zero) + c
        return d

    def mul(c, v):
        return {k: c * x for k, x in v.items()}

    def merge(v, cx):
        d = {}
        for k, c in v.items():
            kk = 'er' if (cx.same and k == 'erp') else k
            d[kk] = d.get(kk, zero) + R(c)
        return d

    def level_of(k, cx):
        return {'er': cx.lr, 'erp': cx.lrp}.get(k)

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
                elif a.startswith('one') or b.startswith('one'):
                    ko, ke = (a, b) if a.startswith('one') else (b, a)
                    i = int(ko[3])
                    val = diagw(i) if level_of(ke, cx) == i else None
                else:
                    val = diagw(level_of(a, cx)) if a == b else None
                if val is not None:
                    tot += c1 * c2 * val
        return tot

    def lin(x, cx, coef):
        tot = zero
        for k, c in x.items():
            if k == 't':
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
        return x.get(key, zero) + x.get(f'one{l}', zero)

    def ipE(x, y, cx, mt, mU, lam, nt, nU):
        x, y = merge(x, cx), merge(y, cx)
        val = x.get('t', zero) * y.get('t', zero) * mt + wdot(x, y, cx, mU)
        nx = x.get('t', zero) * nt + lin(x, cx, nU)
        ny = y.get('t', zero) * nt + lin(y, cx, nU)
        return val - nx * ny * lam

    def cnt(i, j):
        return UU[i] * (UU[j] - (1 if i == j else 0))

    def sum_rr(f):
        tot = zero
        for i in LV:
            tot += f(Ctx(True, i, i)) * UU[i]
            for j in LV:
                if i == j and UU[i] == Kf.one:
                    continue
                tot += f(Ctx(False, i, j)) * cnt(i, j)
        return tot

    def sum_r(f):
        tot = zero
        for i in LV:
            tot += f(Ctx(True, i, i)) * UU[i]
        return tot

    S = sum((UU[i] * SS[i] for i in LV), Kf.zero)
    K = SIG + S / 2
    d = {i: SS[i] + SIG + S / 2 for i in LV}
    lamT = Kf.one / (2 + sum((UU[i] * SS[i] / d[i] for i in LV), Kf.zero))
    d0t = SIG + (SIG + S) / 2
    d0 = {i: SS[i] + (SIG + S) / 2 for i in LV}
    lam0 = Kf.one / (2 + SIG / d0t + sum((UU[i] * SS[i] / d0[i] for i in LV), Kf.zero))

    P = {i: G[f'p{i}'] for i in LV}
    O = {}
    for i in LV:
        for j in LV:
            O[(i, j)] = G[f'o{i}{j}'] if f'o{i}{j}' in G else zero
    q, rho = G['q'], G['rho']
    Z = {i: G[f'z{i}'] for i in LV}

    def Xs(i, j, same):
        return (P[i] if (same and i == j) else O[(i, j)]) * (OM[i] * OM[j])

    def row(i, key):
        v = {}
        for j in LV:
            v = add(v, {f'one{j}': Xs(i, j, False)})
        return add(v, {key: Xs(i, i, True) - Xs(i, i, False)})

    def Tp(psi, key, cx):
        l = level_of('er' if (cx.same and key == 'erp') else key, cx)
        wpsi = lin(merge(psi, cx), cx, lambda i: OM[i])
        c = coord(psi, key, cx)
        e = {key: R.one}
        return add(mul(wpsi * HALF, e), mul(R(OM[l] / 2), psi), mul(c * OM[l], e))

    A = sum_rr(lambda cx: Xs(cx.lr, cx.lrp, cx.same) ** 2) * HALF + sum_r(lambda cx: Xs(cx.lr, cx.lr, True) ** 2)
    trX = sum((Xs(i, i, True) * UU[i] for i in LV), zero)
    wXw = zero
    for i in LV:
        wXw += Xs(i, i, True) * (UU[i] * SS[i])
        for j in LV:
            if i == j and UU[i] == Kf.one:
                continue
            wXw += Xs(i, j, False) * (cnt(i, j) * OM[i] * OM[j])
    zeta = {f'one{i}': Z[i] * (WT * OM[i]) for i in LV}
    wzeta = sum((Z[i] * (UU[i] * OM[i] * WT * OM[i]) for i in LV), zero)
    zz = sum((Z[i] ** 2 * (UU[i] * SIG * SS[i]) for i in LV), zero)
    A += (trX + rho) ** 2 * HALF + rho ** 2 * (3 * HALF) + zz
    cpl = q * 2 * (-trX + wXw * (Kf.one / (2 * K)) - (rho * (2 * K) + wzeta * WT) * (Kf.one / (2 * K - SIG))) + q ** 2 * 2
    wy = rho * WT + wzeta

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
        y = add({'t': rho}, zeta)
        yr = coord(zeta, 'er', cx)
        er = {'er': R.one}
        T1y = add(mul(wy * HALF, er), mul(R(OM[i] / 2), y), mul(yr * OM[i], er))
        return ipE(T2, T1y, cx, Kf.one / d0t, lambda l: Kf.one / d0[l], lam0, WT / d0t, lambda l: OM[l] / d0[l])

    Phi = A + cpl - sum_rr(BU) - sum_r(cross) * 2
    elim = [G[n] for n in 'q p1 p2 p3 o12 o21 o13 o31 o23 o32 o33'.split()]
    keep = [rho, Z[1], Z[2], Z[3]]
    return R, Phi, elim, keep


def to_flint(c, K):
    def conv(p):
        out = {}
        for mon, co in p.terms():
            if mon[0] % 2 or mon[1] % 2 or mon[2] % 2:
                raise ValueError("odd exponent of a square root")
            fq = Fraction(str(co))
            out[(mon[0] // 2, mon[1] // 2, mon[2] // 2, mon[3])] = flint.fmpq(fq.numerator, fq.denominator)
        return K.ctx.from_dict(out)
    return FF(K, conv(c.numer), conv(c.denom))


if __name__ == '__main__':
    t0 = time.time()
    R, Phi, elim, keep = make()
    allv = elim + keep
    idx = {R.gens.index(v): k for k, v in enumerate(allv)}
    n = len(allv)
    Hk = [[Kf.zero] * n for _ in range(n)]
    for mon, co in Phi.terms():
        vs = []
        for pos, e in enumerate(mon):
            vs += [idx[pos]] * e
        assert len(vs) == 2, mon
        a, b = vs
        if a == b:
            Hk[a][a] += co
        else:
            Hk[a][b] += co / 2
            Hk[b][a] += co / 2
    K = Field(['s1', 's2', 's3', 'u3'])
    H = [[to_flint(Hk[i][j], K) for j in range(n)] for i in range(n)]
    print(f'built and converted ({time.time() - t0:.1f} s); eliminating {len(elim)}', flush=True)
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
        print(f'  pivot {p + 1}/{ne} ({time.time() - t0:.1f} s)', flush=True)
    nk = len(keep)
    res = [[H[ne + i][ne + j] for j in range(nk)] for i in range(nk)]
    sym = all(res[i][j] == res[j][i] for i in range(nk) for j in range(nk))
    out = [[(res[i][j].n.to_dict(), res[i][j].d.to_dict()) for j in range(nk)] for i in range(nk)]
    pickle.dump(dict(Q4=out, names=['s1', 's2', 's3', 'u3'], sigma=1), open('three_level.pkl', 'wb'))
    print(f'done ({time.time() - t0:.1f} s); symmetric {sym}; wrote three_level.pkl')
    print('PASS' if sym else 'FAIL')
