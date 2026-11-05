"""Two-level U-weights, every (u1, u2), symbolic (SymPy).  m = 1, T = {t}; U = U1 + U2 with |U_i| = u_i, w_h = om_i on
U_i, s_i = om_i^2, sig = s_t.  Gauge x_t = 0; the eliminated block is the (C_u) block, positive definite (Theorem C),
and invariant under S_{u1} x S_{u2}, so for invariant letter-2 data (rho, z_h = z_i on U_i) the unique minimiser is
invariant:  x_rr = p_i (r in U_i),  x_rh = o_ij (r in U_i, h in U_j, r != h)  (scaled: X_rh = w_r w_h x_rh), and q.
The closed form is the one of equal_u.py (cf.py with x_t = 0).

Cases: 'gen' (u1, u2 >= 2 symbolic), 'one' (u1 = 1, u2 >= 2 symbolic; o_11 absent).
Output: two_level_{case}.pkl with the 3 x 3 matrix Q3 of G_z on the invariant letter-2 data (rho, z1, z2):
  Q3[0,0] = G_z[rho,rho], Q3[0,i] = sum_{h in U_i} G_z[rho,h] = u_i G_z[rho,h_i],
  Q3[i,j] = sum_{h in U_i, h' in U_j} G_z[h,h'].
"""
import sys, pickle, time
import sympy as sp
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from flintfield import Field, F as FF

sig = sp.symbols('sigma', positive=True)
s1, s2 = sp.symbols('s1 s2', positive=True)
u1, u2 = sp.symbols('u1 u2', positive=True)
wt, o1, o2 = sp.sqrt(sig), sp.sqrt(s1), sp.sqrt(s2)
OM = {1: o1, 2: o2}
SS = {1: s1, 2: s2}
UU = {1: u1, 2: u2}


class Ctx:
    def __init__(self, same, lr, lrp):
        self.same, self.lr, self.lrp = same, lr, (lr if same else lrp)


def V(**kw):
    return {k: sp.sympify(v) for k, v in kw.items()}


def add(*vs):
    d = {}
    for v in vs:
        for k, c in v.items():
            d[k] = d.get(k, 0) + c
    return d


def mul(c, v):
    return {k: c * x for k, x in v.items()}


def merge(v, cx):
    if not cx.same:
        return v
    d = {}
    for k, c in v.items():
        kk = 'er' if k == 'erp' else k
        d[kk] = d.get(kk, 0) + c
    return d


def level_of(k, cx):
    return {'er': cx.lr, 'erp': cx.lrp}.get(k)


def wdot(x, y, cx, diagw):
    """sum_h x_h y_h diagw(level h) over U (x, y merged)."""
    tot = 0
    for k1, c1 in x.items():
        if k1 == 't':
            continue
        for k2, c2 in y.items():
            if k2 == 't':
                continue
            if k1.startswith('one') and k2.startswith('one'):
                i, j = int(k1[3]), int(k2[3])
                val = UU[i] * diagw(i) if i == j else 0
            elif k1.startswith('one') or k2.startswith('one'):
                ko, ke = (k1, k2) if k1.startswith('one') else (k2, k1)
                i = int(ko[3])
                val = diagw(i) if level_of(ke, cx) == i else 0
            else:
                val = diagw(level_of(k1, cx)) if k1 == k2 else 0
            tot += c1 * c2 * val
    return tot


def lin(x, cx, coef):
    """sum_h coef(level h) x_h over U."""
    tot = 0
    for k, c in x.items():
        if k == 't':
            continue
        if k.startswith('one'):
            i = int(k[3])
            tot += c * UU[i] * coef(i)
        else:
            tot += c * coef(level_of(k, cx))
    return tot


def coord(x, key, cx):
    """coordinate of x at the copy named by key ('er' or 'erp')."""
    x = merge(x, cx)
    key = 'er' if (cx.same and key == 'erp') else key
    l = level_of(key, cx)
    return x.get(key, 0) + x.get(f'one{l}', 0)


def ipE(x, y, cx, mt, mU, lam, nt, nU):
    """<x, M y>, M = diag(mt on t, mU(level) on U) - lam nu nu^T, nu = nt e_t + nU(level) on U."""
    x, y = merge(x, cx), merge(y, cx)
    val = mt * x.get('t', 0) * y.get('t', 0) + wdot(x, y, cx, mU)
    nx = nt * x.get('t', 0) + lin(x, cx, nU)
    ny = nt * y.get('t', 0) + lin(y, cx, nU)
    return val - lam * nx * ny


def build(case):
    p1, p2, o11, o12, o21, o22, q, rho, z1, z2 = sp.symbols('p1 p2 o11 o12 o21 o22 q rho z1 z2')
    P = {1: p1, 2: p2}
    O = {(1, 1): o11, (1, 2): o12, (2, 1): o21, (2, 2): o22}
    if case == 'one':
        UU[1] = sp.Integer(1)
        O[(1, 1)] = sp.Integer(0)
    else:
        UU[1] = u1
    uA, uB = UU[1], UU[2]
    S = uA * s1 + uB * s2
    K = sig + S / 2

    def Xs(i, j, same):   # unscaled entry X_rh, r in U_i, h in U_j
        return OM[i] * OM[j] * (P[i] if (same and i == j) else O[(i, j)])

    def row(i, key):
        """row r (r in U_i) as a vector on U: off-diagonal value on every copy, diagonal correction at r."""
        v = add(mul(Xs(i, 1, False), V(one1=1)), mul(Xs(i, 2, False), V(one2=1)))
        return add(v, {key: Xs(i, i, True) - Xs(i, i, False)})

    def sum_rr(f):
        """sum over r, r' in U of f(cx) (f returns an expression)."""
        tot = 0
        for i in (1, 2):
            tot += UU[i] * f(Ctx(True, i, i))
            for j in (1, 2):
                cnt = UU[i] * (UU[j] - (1 if i == j else 0))
                if cnt != 0:
                    tot += cnt * f(Ctx(False, i, j))
        return tot

    def sum_r(f):
        return sum((UU[i] * f(Ctx(True, i, i)) for i in (1, 2)), sp.Integer(0))

    # A-part
    A = sp.Rational(1, 2) * sum_rr(lambda cx: (Xs(cx.lr, cx.lrp, cx.same)) ** 2) + sum_r(lambda cx: Xs(cx.lr, cx.lr, True) ** 2)
    trX = sum((UU[i] * Xs(i, i, True) for i in (1, 2)), sp.Integer(0))
    # w^T X w
    wXw = 0
    for i in (1, 2):
        wXw += UU[i] * OM[i] ** 2 * Xs(i, i, True)
        for j in (1, 2):
            wXw += UU[i] * (UU[j] - (1 if i == j else 0)) * OM[i] * OM[j] * Xs(i, j, False)
    zeta = add(mul(wt * o1 * z1, V(one1=1)), mul(wt * o2 * z2, V(one2=1)))
    wzeta = UU[1] * o1 * wt * o1 * z1 + UU[2] * o2 * wt * o2 * z2
    zz = UU[1] * (wt * o1 * z1) ** 2 + UU[2] * (wt * o2 * z2) ** 2
    A += sp.Rational(1, 2) * (trX + rho) ** 2 + sp.Rational(3, 2) * rho ** 2 + zz
    cpl = 2 * q * (-trX + wXw / (2 * K) - (2 * K * rho + wt * wzeta) / (2 * K - sig)) + 2 * q ** 2
    # B_U with E_T = diag(d_i) + 1/2 w w^T
    d = {i: SS[i] + sig + S / 2 for i in (1, 2)}
    gT = sum((UU[i] * SS[i] / d[i] for i in (1, 2)), sp.Integer(0))
    lamT = 1 / (2 + gT)

    def Tp(psi, key, cx):
        """T_{r_key} psi = 1/2 (w.psi) e + 1/2 w_key psi + w_key psi_key e."""
        l = level_of('er' if (cx.same and key == 'erp') else key, cx)
        wpsi = lin(merge(psi, cx), cx, lambda i: OM[i])
        c = coord(psi, key, cx)
        e = {key: sp.Integer(1)}
        return add(mul(sp.Rational(1, 2) * wpsi, e), mul(OM[l] / 2, psi), mul(OM[l] * c, e))

    def BU(cx):
        psi_r = row(cx.lr, 'er')
        psi_rp = row(cx.lrp, 'erp') if not cx.same else row(cx.lr, 'er')
        a1 = Tp(psi_r, 'erp', cx)
        a2 = Tp(psi_rp, 'er', cx)
        return ipE(a1, a2, cx, 0, lambda i: 1 / d[i], lamT, 0, lambda i: OM[i] / d[i])
    BUv = sum_rr(BU)
    # cross term with E_0 on U + t
    d0t = sig + (sig + S) / 2
    d0 = {i: SS[i] + (sig + S) / 2 for i in (1, 2)}
    g0 = sig / d0t + sum((UU[i] * SS[i] / d0[i] for i in (1, 2)), sp.Integer(0))
    lam0 = 1 / (2 + g0)

    def cross(cx):
        i = cx.lr
        psi_r = row(i, 'er')
        wpsi = lin(psi_r, cx, lambda l: OM[l])
        T2 = add(mul(sp.Rational(1, 2) * wpsi, V(t=1)), mul(-wt / 2, psi_r))
        y = add(V(t=rho), zeta)
        wy = wt * rho + wzeta
        yr = coord(zeta, 'er', cx)
        er = V(er=1)
        T1y = add(mul(sp.Rational(1, 2) * wy, er), mul(OM[i] / 2, y), mul(OM[i] * yr, er))
        return ipE(T2, T1y, cx, 1 / d0t, lambda l: 1 / d0[l], lam0, wt / d0t, lambda l: OM[l] / d0[l])
    crossv = sum_r(cross)
    Phi = sp.expand(A + cpl - BUv - 2 * crossv)
    elim = [p1, p2, o12, o21, o22, q] + ([] if case == 'one' else [o11])
    keep = [rho, z1, z2]
    return Phi, elim, keep


def to_flint(e, K, gens):
    num, den = sp.fraction(sp.together(e))
    pn = sp.Poly(sp.expand(num), *gens)
    pd = sp.Poly(sp.expand(den), *gens)
    cn = K.ctx.from_dict({m: flint.fmpq(int(sp.Rational(c).p), int(sp.Rational(c).q)) for m, c in pn.as_dict().items()})
    cd = K.ctx.from_dict({m: flint.fmpq(int(sp.Rational(c).p), int(sp.Rational(c).q)) for m, c in pd.as_dict().items()})
    return FF(K, cn, cd)


def derive(case):
    t0 = time.time()
    Phi, elim, keep = build(case)
    for pw in Phi.atoms(sp.Pow):
        if pw.base in (sig, s1, s2) and not pw.exp.is_integer:
            raise ValueError(f'non-integer exponent {pw}')
    allv = elim + keep
    gens = [sig, s1, s2, u2] + ([] if case == 'one' else [u1])
    names = [str(g) for g in gens]
    K = Field(names)
    n = len(allv)
    H = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            e = sp.diff(Phi, allv[i], allv[j]) / 2
            H[i][j] = H[j][i] = to_flint(e, K, gens)
    print(f'{case}: form built and converted ({time.time() - t0:.1f} s), {len(elim)} pivots', flush=True)
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
    Q3 = [[H[ne + i][ne + j] for j in range(3)] for i in range(3)]
    return Q3, names, Phi


if __name__ == '__main__':
    case = sys.argv[1]
    Q3, names, Phi = derive(case)
    out = [[(Q3[i][j].n.to_dict(), Q3[i][j].d.to_dict()) for j in range(3)] for i in range(3)]
    sym = all(Q3[i][j] == Q3[j][i] for i in range(3) for j in range(3))
    pickle.dump(dict(Q3=out, names=names, case=case), open(f'two_level_{case}_exprtree.pkl', 'wb'))
    print('symmetric:', sym, '; wrote', f'two_level_{case}_exprtree.pkl')
    print('PASS' if sym else 'FAIL')
