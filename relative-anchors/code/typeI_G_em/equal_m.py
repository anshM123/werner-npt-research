"""Equal U-weights, every u >= 3 and every number m of copies carrying the letter 2: the blocks G_tt and G_tt' of the
Schur complement G, symbolically in u (SymPy for the forms, python-flint for the eliminations).

Setting (relative-anchors/README.md Section 3; typeI_G_u2/NOTES.md Section 1).  U = u copies of equal weight s = om^2,
T = copies carrying one letter 2, sigma = s_T, S = u s, K = sigma + S/2.  Work in the gauge sum_t w_t psi_t = 0
(psi_t = row t of X).  Then (typeI_G_u2/NOTES.md, Section 1)
    G_tt  = diag(2, I) - (1/4)(1 - s_t/sigma) L_t^T L_t - C_t^T C^{-1} C_t,
    G_tt' = (1/2) e_rho e_rho^T - P_tt' + (w_t w_t'/(4 sigma)) L_t^T L_t' - C_t^T C^{-1} C_t'      (t != t'),
where C is the (C_u) block on (X_U, q) (Theorem C), C_t the coupling of y_t = (rho_t, zeta_t) to (X_U, q), L_t the map
y_t -> l_t = -zeta_t - 2 T_(t,2)^T E_{T-t}^{-1} (T_(t,1) - w_t) y_t, and P_tt' the pair term
<T_(t',2) y_t, E_{T-tt'}^{-1} T_(t,2) y_t'>.  Only (s_t, s_t', sigma, s, u) enter.

By the permutations of U (C positive definite, so the minimiser is unique and equivariant) the blocks split into
  * the invariant part: X_U = s (a I + b (J - I)), q, zeta_t = w_t om z_t 1_U   (coordinates (rho_t, z_t) per copy);
  * the standard part:  X_U = s (al diag f + be f 1^T + ga 1 f^T), q = 0, rho = 0, zeta_t = w_t om z_t f, sum f = 0
    (u >= 3: the standard representation occurs three times in R^{U x U}); value per |f|^2.
Here the configuration T = {t1, t2} (+ further copies of zero data) is used with symbols s1, s2 and total sigma >= s1 + s2.
Output (pickled, as strings of flint polynomials): Ginv (4 x 4 on (rho1, z1, rho2, z2)) and Gstd (2 x 2 on (z1, z2)),
entries rational in (u, sigma, s, s1, s2); G_tt is the (t, t) block, G_tt' the (t, t') block.
"""
import sys, os, pickle, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
import equal_u as EU
from equal_u import Vec, dotU, sumU, ip, sum_r, u, sig, s

s1, s2 = sp.symbols('s1 s2', positive=True)
w1, w2, om = sp.sqrt(s1), sp.sqrt(s2), sp.sqrt(s)
S = u * s
K = sig + S / 2
half = sp.Rational(1, 2)


def coupling(row, zeta, rho, q, trX, wi, si):
    """bilinear terms between (X_U, q) and y_i = rho e_t + zeta for one copy t of weight si = wi^2:
       trX rho  - 2 q (2K rho + wi w_U.zeta)/(2K - si)  - 2 sum_{r in U} <T_(t,2) psi_r, E_{T-t}^{-1} T_(r,1) y>."""
    dt = K + si / 2
    dU = s + K - si / 2
    g0 = si / dt + u * s / dU
    lam0 = 1 / (2 + g0)
    nt, nU = wi / dt, om / dU
    y = Vec({'t': rho}) + zeta
    wy = wi * rho + om * sumU(zeta)
    psi_r = row(True)[0]
    er = Vec({'er': 1})
    T2 = (half * om * sumU(psi_r)) * Vec({'t': 1}) + (-wi / 2) * psi_r
    yr = dotU(zeta, er, True)
    T1y = (half * wy) * er + (om / 2) * y + (om * yr) * er
    cross = ip(T2, T1y, 1 / dt, 1 / dU, lam0, nt, nU, True)
    wzeta = om * sumU(zeta)
    return trX * rho - 2 * q * (2 * K * rho + wi * wzeta) / (2 * K - si) - 2 * sum_r(cross)


def ell_inv(rho, c, wi, si):
    """coefficient of 1_U in l_i, for y = rho e_t + c 1_U (invariant data)."""
    dt = K + si / 2
    dU = s + K - si / 2
    g = si / dt + u * s / dU
    vt = wi * rho + half * om * u * c
    vU = -half * wi * c
    ev = (wi / dt) * vt + u * (om / dU) * vU
    Et = vt / dt - (wi / dt) * ev / (2 + g)
    EU_ = vU / dU - (om / dU) * ev / (2 + g)
    return -c - Et * om + wi * EU_


def ell_std(c, wi, si):
    """coefficient of f in l_i for y = c f (standard data)."""
    dU = s + K - si / 2
    return -c * (1 + si / (2 * dU))


def pair_inv(r1, c1, r2, c2):
    """<T_(2,2) y_1, E_{T-12}^{-1} T_(1,2) y_2> for y_i = r_i e_i + c_i 1_U."""
    sh = (s1 + s2) / 2
    d1, d2, dU = s1 + K - sh, s2 + K - sh, s + K - sh
    g = s1 / d1 + s2 / d2 + u * s / dU
    A1, A2, AU = -half * w2 * r1, half * (w1 * r1 + om * u * c1), -half * w2 * c1
    B1, B2, BU = half * (w2 * r2 + om * u * c2), -half * w1 * r2, -half * w1 * c2
    eA = (w1 / d1) * A1 + (w2 / d2) * A2 + u * (om / dU) * AU
    eB = (w1 / d1) * B1 + (w2 / d2) * B2 + u * (om / dU) * BU
    return A1 * B1 / d1 + A2 * B2 / d2 + u * AU * BU / dU - eA * eB / (2 + g)


def pair_std(c1, c2):
    sh = (s1 + s2) / 2
    dU = s + K - sh
    return w1 * w2 * c1 * c2 / (4 * dU)


def derive():
    a, b, q = sp.symbols('a b q')
    al, be, ga = sp.symbols('alpha beta gamma')
    r1, z1, r2, z2 = sp.symbols('rho1 z1 rho2 z2')
    out = {}
    # ---------------- invariant part
    def row_inv(same):
        return Vec({'one': s * b, 'er': s * (a - b)}), Vec({'one': s * b, 'erp': s * (a - b)})

    def Xrr_inv(same):
        return (s * a, s * a) if same else (s * b, s * b)
    trX = u * s * a
    PhiU = EU.form(row_inv, Xrr_inv, Vec(), 0, q, trX=trX, wXw=s * s * u * (a + (u - 1) * b), wzeta=0, zz=0)
    c1, c2 = w1 * om * z1, w2 * om * z2
    cp1 = coupling(row_inv, Vec({'one': c1}), r1, q, trX, w1, s1)
    cp2 = coupling(row_inv, Vec({'one': c2}), r2, q, trX, w2, s2)
    l1, l2 = ell_inv(r1, c1, w1, s1), ell_inv(r2, c2, w2, s2)
    yy = 2 * r1 ** 2 + 2 * r2 ** 2 + r1 * r2 + u * (c1 ** 2 + c2 ** 2)
    yy += -sp.Rational(1, 4) * u * (l1 ** 2 + l2 ** 2) + u * (w1 * l1 + w2 * l2) ** 2 / (4 * sig)
    yy += -2 * pair_inv(r1, c1, r2, c2)
    Phi = PhiU + cp1 + cp2 + yy
    out['inv'] = (Phi, [a, b, q], [r1, z1, r2, z2])
    # ---------------- standard part (per |f|^2)
    F2 = EU.F2
    fr, frp = EU.fr, EU.frp

    def row_std(same):
        pr = Vec({'er': s * al * fr, 'one': s * be * fr, 'f': s * ga})
        prp = Vec({'erp': s * al * frp, 'one': s * be * frp, 'f': s * ga})
        return pr, prp

    def Xrr_std(same):
        if same:
            return (s * (al * fr + be * fr + ga * fr),) * 2
        return (s * (be * fr + ga * frp), s * (be * frp + ga * fr))
    PhiUs = EU.form(row_std, Xrr_std, Vec(), 0, 0, trX=0, wXw=0, wzeta=0, zz=0)
    cs1 = coupling(row_std, Vec({'f': c1}), 0, 0, 0, w1, s1)
    cs2 = coupling(row_std, Vec({'f': c2}), 0, 0, 0, w2, s2)
    ls1, ls2 = ell_std(c1, w1, s1), ell_std(c2, w2, s2)
    yys = (c1 ** 2 + c2 ** 2) * F2 - sp.Rational(1, 4) * (ls1 ** 2 + ls2 ** 2) * F2 \
        + (w1 * ls1 + w2 * ls2) ** 2 * F2 / (4 * sig) - 2 * pair_std(c1, c2) * F2
    Phis = sp.expand((PhiUs + cs1 + cs2 + yys) / F2)
    assert not Phis.has(F2), 'standard form not proportional to |f|^2'
    out['std'] = (Phis, [al, be, ga], [z1, z2])
    return out


def to_flint(Phi, elim, keep, names=('u', 'sigma', 's', 's1', 's2')):
    """Hessian/2 of the quadratic form (second derivatives, no expansion), converted exactly into the flint rational
    function field: each square root sqrt(x) of a weight becomes an auxiliary generator W_x, and at the end every
    exponent of every W_x must be even (else ValueError), W_x^(2k) -> x^k.  Then the Schur complement onto `keep`."""
    from flintfield import Field
    wnames = ['W_' + n for n in names]
    Fd8 = Field(list(names) + wnames)
    g8 = dict(zip(list(names) + wnames, Fd8.gens))
    Fd = Field(list(names))
    nv = len(names)

    def conv(e):
        if e.is_Rational:
            return Fd8(int(e.p)) / Fd8(int(e.q))
        if e.is_Symbol:
            return g8[e.name]
        if e.is_Add:
            r = Fd8(0)
            for a in e.args:
                r = r + conv(a)
            return r
        if e.is_Mul:
            r = Fd8(1)
            for a in e.args:
                r = r * conv(a)
            return r
        if e.is_Pow:
            b, x = e.args
            if x.is_Integer:
                base = conv(b)
                n = int(x)
                r = Fd8(1)
                for _ in range(abs(n)):
                    r = r * base
                return r if n >= 0 else Fd8(1) / r
            if x.is_Rational and x.q == 2 and b.is_Symbol:
                n = int(x.p)
                base = g8['W_' + b.name]
                r = Fd8(1)
                for _ in range(abs(n)):
                    r = r * base
                return r if n >= 0 else Fd8(1) / r
        raise ValueError('cannot convert %s' % e)

    def deW(p):
        d = {}
        for mon, c in zip(p.monoms(), p.coeffs()):
            mon = list(mon)
            for i in range(nv):
                if mon[nv + i] % 2:
                    raise ValueError('odd power of a weight left')
                mon[i] += mon[nv + i] // 2
            d[tuple(mon[:nv])] = d.get(tuple(mon[:nv]), 0) + c
        return Fd.ctx.from_dict(d)

    from flintfield import F as FE
    allv = elim + keep
    n = len(allv)
    Hf = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            c = sp.diff(Phi, allv[i], allv[j]) / 2
            c = c.subs({v: 0 for v in allv})
            x8 = conv(c)
            Hf[i][j] = Hf[j][i] = FE(Fd, deW(x8.n), deW(x8.d))
    ne = len(elim)
    A = [row[:] for row in Hf]
    for p in range(ne):
        pv = A[p][p]
        for i in range(p + 1, n):
            if A[i][p] == 0:
                continue
            f = A[i][p] / pv
            for j in range(p + 1, n):
                A[i][j] = A[i][j] - f * A[p][j]
    G = [[A[ne + i][ne + j] for j in range(len(keep))] for i in range(len(keep))]
    return Fd, Hf, G


def ser(x):
    """exact serialisation of a flint rational function: (num, den) as dicts monomial -> (p, q)."""
    f = lambda p: {tuple(int(e) for e in m): (int(c.p), int(c.q)) for m, c in zip(p.monoms(), p.coeffs())}
    return (f(x.n), f(x.d))


if __name__ == '__main__':
    t0 = time.time()
    R = derive()
    print('forms built in %.1f s' % (time.time() - t0), flush=True)
    res = {}
    for key in ('inv', 'std'):
        Phi, elim, keep = R[key]
        Fd, Hf, G = to_flint(Phi, elim, keep)
        res[key] = dict(G=[[ser(x) for x in row] for row in G], H=[[ser(x) for x in row] for row in Hf],
                        keep=[str(v) for v in keep], elim=[str(v) for v in elim])
        print(key, 'eliminated in %.1f s' % (time.time() - t0), flush=True)
        for i in range(len(keep)):
            for j in range(i, len(keep)):
                print('  G[%s,%s]: num terms %d, den terms %d' % (keep[i], keep[j], len(G[i][j].n), len(G[i][j].d)))
    pickle.dump(res, open('equal_m.pkl', 'wb'))
    print('done in %.1f s' % (time.time() - t0))
    print('PASS')
