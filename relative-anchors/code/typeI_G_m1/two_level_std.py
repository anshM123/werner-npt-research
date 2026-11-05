"""Two-level U-weights: the standard part of level 1 (SymPy, symbolic in the multiplicities).
U = U1 + U2, |U_i| = u_i, w_h = om_i on U_i, sig = s_t.  Gauge x_t = 0.  For letter-2 data zeta = w_t om_1 z f with f a
vector on U1, sum f = 0 (rho = 0), the unique minimiser of the (C_u) block is S_{u1} x S_{u2}-equivariant in f, hence
(for u1 >= 3) of the form (scaled coordinates X_rh = w_r w_h x_rh)
  x_rr = (al + be + ga) f_r,  x_rh = be f_r + ga f_h  (r != h in U1),
  x_rh = de f_r (r in U1, h in U2),  x_hr = ep f_r (h in U2, r in U1),  x = 0 on U2 x U2,  q = 0.
(The standard representation of S_{u1} occurs 3 times in the U1 x U1 block for u1 >= 3, once in each of the U1 x U2
and U2 x U1 blocks, and not in the U2 x U2 block or in q.)  The Schur complement onto z is lam1 |f|^2 z^2 with
  lam1 = G_z[h,h] - G_z[h,h'],  h != h' in U1.
Cases: 'gen' (u1 >= 3, u2 >= 1 symbolic), 'eq' is not needed (equal_u.py).  Output two_level_std.pkl."""
import sys, os, pickle, time
import sympy as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from flintfield import Field, F as FF

sig = sp.symbols('sigma', positive=True)
s1, s2 = sp.symbols('s1 s2', positive=True)
u1, u2 = sp.symbols('u1 u2', positive=True)
wt, o1, o2 = sp.sqrt(sig), sp.sqrt(s1), sp.sqrt(s2)
fr, frp, F2 = sp.symbols('f_r f_rp F2')
OM = {1: o1, 2: o2}
SS = {1: s1, 2: s2}
UU = {1: u1, 2: u2}


class Ctx:
    def __init__(self, same, lr, lrp):
        self.same, self.lr, self.lrp = same, lr, (lr if same else lrp)


def add(*vs):
    d = {}
    for v in vs:
        for k, c in v.items():
            d[k] = d.get(k, 0) + c
    return d


def mul(c, v):
    return {k: c * x for k, x in v.items()}


def fix(e, cx):
    """f vanishes on U2: substitute f_r = 0 (r in U2), f_r' = 0 (r' in U2); identify r' = r if same."""
    e = sp.sympify(e)
    if cx.same:
        e = e.subs(frp, fr)
    if cx.lr == 2:
        e = e.subs(fr, 0)
    if cx.lrp == 2:
        e = e.subs(frp, 0)
    return e


def merge(v, cx):
    d = {}
    for k, c in v.items():
        kk = 'er' if (cx.same and k == 'erp') else k
        d[kk] = d.get(kk, 0) + fix(c, cx)
    return d


def level_of(k, cx):
    return {'er': cx.lr, 'erp': cx.lrp}.get(k)


def fval(k, cx):
    """value of f at the copy named k"""
    return {'er': fix(fr, cx), 'erp': fix(frp, cx)}[k]


def wdot(x, y, cx, diagw):
    tot = 0
    for k1, c1 in x.items():
        if k1 == 't':
            continue
        for k2, c2 in y.items():
            if k2 == 't':
                continue
            a, b = sorted((k1, k2))
            if a.startswith('one') and b.startswith('one'):
                i, j = int(a[3]), int(b[3])
                val = UU[i] * diagw(i) if i == j else 0
            elif a == 'f' and b == 'f':
                val = F2 * diagw(1)
            elif (a.startswith('one') and b == 'f') or (b.startswith('one') and a == 'f'):
                val = 0
            elif a.startswith('one') or b.startswith('one'):
                ko, ke = (a, b) if a.startswith('one') else (b, a)
                i = int(ko[3])
                val = diagw(i) if level_of(ke, cx) == i else 0
            elif a == 'f' or b == 'f':
                ke = b if a == 'f' else a
                val = fval(ke, cx) * diagw(1)
            else:
                val = diagw(level_of(a, cx)) if a == b else 0
            tot += c1 * c2 * val
    return tot


def lin(x, cx, coef):
    tot = 0
    for k, c in x.items():
        if k == 't' or k == 'f':
            continue          # sum of f is 0
        if k.startswith('one'):
            i = int(k[3])
            tot += c * UU[i] * coef(i)
        else:
            tot += c * coef(level_of(k, cx))
    return tot


def coord(x, key, cx):
    x = merge(x, cx)
    key = 'er' if (cx.same and key == 'erp') else key
    l = level_of(key, cx)
    return x.get(key, 0) + x.get(f'one{l}', 0) + x.get('f', 0) * fval(key, cx)


def ipE(x, y, cx, mt, mU, lam, nt, nU):
    x, y = merge(x, cx), merge(y, cx)
    val = mt * x.get('t', 0) * y.get('t', 0) + wdot(x, y, cx, mU)
    nx = nt * x.get('t', 0) + lin(x, cx, nU)
    ny = nt * y.get('t', 0) + lin(y, cx, nU)
    return val - lam * nx * ny


def red_pairs(e, cx):
    """sum over the copies of the context of e (polynomial of degree <= 2 in f_r, f_r', F2, after fix)."""
    e = sp.expand(fix(e, cx))
    p = sp.Poly(e, fr, frp, F2)
    i, j = cx.lr, cx.lrp
    tot = 0
    for (a, b, c), co in p.terms():
        deg = a + b + 2 * c
        if deg not in (0, 2):
            if deg == 1:
                continue       # sums of f vanish
            raise ValueError('degree > 2 in f')
        if cx.same:
            # sum over r in U_i
            if c == 1:
                val = UU[i] * F2
            elif (a, b) == (0, 0):
                val = UU[i]
            elif a + b == 2:
                val = F2 if i == 1 else 0
            else:
                val = 0
        else:
            cnt = UU[i] * (UU[j] - (1 if i == j else 0))
            if c == 1:
                val = cnt * F2
            elif (a, b) == (0, 0):
                val = cnt
            elif (a, b) == (2, 0):      # sum_{r in U_i, r' in U_j, r != r'} f_r^2
                val = (UU[j] - (1 if i == j else 0)) * F2 if i == 1 else 0
            elif (a, b) == (0, 2):
                val = (UU[i] - (1 if i == j else 0)) * F2 if j == 1 else 0
            elif (a, b) == (1, 1):
                val = -F2 if (i == 1 and j == 1) else 0
            else:
                val = 0
        tot += co * val
    return tot


def build(case):
    al, be, ga, de, ep, z = sp.symbols('alpha beta gamma delta epsilon z')
    if case == 'u1is2':
        UU[1] = sp.Integer(2)
        ga = sp.Integer(0)
    else:
        UU[1] = u1
    S = UU[1] * s1 + UU[2] * s2
    K = sig + S / 2

    def xs(i, j, same, a, b):
        """unscaled X_rh for r in U_i (f-value a = f_r), h in U_j (f-value b = f_h)."""
        pre = OM[i] * OM[j]
        if i == 1 and j == 1:
            return pre * ((al + be + ga) * a if same else be * a + ga * b)
        if i == 1 and j == 2:
            return pre * de * a
        if i == 2 and j == 1:
            return pre * ep * b
        return 0

    def row(i, key):
        """row r (r in U_i) as a vector on U."""
        a = fr if key == 'er' else frp
        if i == 1:
            v = {'one1': OM[1] * OM[1] * be * a, 'f': OM[1] * OM[1] * ga, 'one2': OM[1] * OM[2] * de * a,
                 key: OM[1] * OM[1] * al * a}
        else:
            v = {'f': OM[2] * OM[1] * ep}
        return v

    def sum_rr(f):
        tot = 0
        for i in (1, 2):
            cx = Ctx(True, i, i)
            tot += red_pairs(f(cx), cx)
            for j in (1, 2):
                cx = Ctx(False, i, j)
                tot += red_pairs(f(cx), cx)
        return tot

    def sum_r(f):
        return sum((red_pairs(f(Ctx(True, i, i)), Ctx(True, i, i)) for i in (1, 2)), sp.Integer(0))

    A = sp.Rational(1, 2) * sum_rr(lambda cx: xs(cx.lr, cx.lrp, cx.same, fr, frp if not cx.same else fr) ** 2)
    A += sum_r(lambda cx: xs(cx.lr, cx.lr, True, fr, fr) ** 2)
    zeta = {'f': wt * o1 * z}
    A += (wt * o1 * z) ** 2 * F2
    # trX = 0, w^T X w = 0, w.zeta = 0 for this data, rho = q = 0: the q-row and the trace terms vanish
    d = {i: SS[i] + sig + S / 2 for i in (1, 2)}
    gT = sum((UU[i] * SS[i] / d[i] for i in (1, 2)), sp.Integer(0))
    lamT = 1 / (2 + gT)

    def Tp(psi, key, cx):
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
    d0t = sig + (sig + S) / 2
    d0 = {i: SS[i] + (sig + S) / 2 for i in (1, 2)}
    g0 = sig / d0t + sum((UU[i] * SS[i] / d0[i] for i in (1, 2)), sp.Integer(0))
    lam0 = 1 / (2 + g0)

    def cross(cx):
        i = cx.lr
        psi_r = row(i, 'er')
        wpsi = lin(merge(psi_r, cx), cx, lambda l: OM[l])
        T2 = add(mul(sp.Rational(1, 2) * wpsi, {'t': 1}), mul(-wt / 2, psi_r))
        y = dict(zeta)
        wy = 0                      # w_t rho + w.zeta = 0
        yr = coord(zeta, 'er', cx)
        er = {'er': sp.Integer(1)}
        T1y = add(mul(sp.Rational(1, 2) * wy, er), mul(OM[i] / 2, y), mul(OM[i] * yr, er))
        return ipE(T2, T1y, cx, 1 / d0t, lambda l: 1 / d0[l], lam0, wt / d0t, lambda l: OM[l] / d0[l])
    crossv = sum_r(cross)
    Phi = A - BUv - 2 * crossv
    return Phi, ([al, be, de, ep] if case == 'u1is2' else [al, be, ga, de, ep]), [z]


def to_flint(e, K, gens):
    num, den = sp.fraction(sp.together(e))
    pn = sp.Poly(sp.expand(num), *gens)
    pd = sp.Poly(sp.expand(den), *gens)
    cv = lambda c: flint.fmpq(int(sp.Rational(c).p), int(sp.Rational(c).q))
    return FF(K, K.ctx.from_dict({m: cv(c) for m, c in pn.as_dict().items()}),
              K.ctx.from_dict({m: cv(c) for m, c in pd.as_dict().items()}))


if __name__ == '__main__':
    t0 = time.time()
    case = sys.argv[1]
    Phi, elim, keep = build(case)
    allv = elim + keep
    gens = [sig, s1, s2, u2] + ([] if case == 'u1is2' else [u1])
    K = Field([str(g) for g in gens])
    n = len(allv)
    H = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            e = sp.diff(Phi, allv[i], allv[j]) / 2
            e = sp.expand(e)
            # must be proportional to F2, with no square roots left
            e = sp.expand(e / F2)
            assert not e.has(F2) and not e.has(fr) and not e.has(frp), 'not proportional to |f|^2'
            for pw in e.atoms(sp.Pow):
                if pw.base in (sig, s1, s2) and not pw.exp.is_integer:
                    raise ValueError(f'non-integer exponent {pw}')
            H[i][j] = H[j][i] = to_flint(e, K, gens)
    print(f'built ({time.time() - t0:.1f} s)', flush=True)
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
    lam1 = H[n - 1][n - 1]
    pickle.dump(dict(lam1=(lam1.n.to_dict(), lam1.d.to_dict()), names=[str(g) for g in gens]),
                open(f'two_level_std_{case}_exprtree.pkl', 'wb'))
    print(f'done ({time.time() - t0:.1f} s); wrote two_level_std_{case}_exprtree.pkl')
    print('PASS')
