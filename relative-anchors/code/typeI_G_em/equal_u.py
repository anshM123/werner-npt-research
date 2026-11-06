"""Equal U-weights, every u, symbolic in u (SymPy).  m = 1, T = {t}, U = {1, ..., u}, w_h = om for h in U,
s = om^2, sig = s_t = w_t^2.  Gauge x_t = 0 (the deleted coordinates P_{t,h}); the eliminated block is the (C_u) block
on (X_U, q), positive definite by Theorem C, and invariant under the permutations of U (equal U-weights), so the unique
minimiser for S_u-invariant (resp. standard) letter-2 data lies in the invariant (resp. matching standard) part.

The closed form of the type I block (cf.py; README Lemma B, N; typeI_u2/NOTES.md Section A1) with x_t = 0 reads
  Phi = 1/2 sum_{r,h} X_rh^2 + sum_h X_hh^2 + 1/2 (tr X + rho)^2 + 3/2 rho^2 + |zeta|^2
        + 2 q [ -tr X + w^T X w/(2K) - (2K rho + w_t w.zeta)/(2K - sig) ] + 2 q^2
        - sum_{r,r' in U} < T_{r'} psi_r, E_T^{-1} T_r psi_{r'} >
        - 2 sum_{r in U} < T_{(t,2)} psi_r, E_0^{-1} T_{(r,1)} y >,
  K = sig + S/2, S = u s, y = rho e_t + zeta (on U + t), psi_r = row r of X,
  T_{r'} psi = 1/2 (w.psi) e_{r'} + 1/2 w_{r'} psi + w_{r'} psi_{r'} e_{r'}  (on U),
  T_{(t,2)} psi = 1/2 (w.psi) e_t - 1/2 w_t psi  (on U + t),
  T_{(r,1)} y = 1/2 (wt.y) e_r + 1/2 w_r y + w_r y_r e_r  (on U + t, wt = (w_t, w_U)),
  E_T = diag(s_h + sig + S/2) + 1/2 w w^T on U;  E_0 = diag(s_x + (sig + S)/2) + 1/2 wt wt^T on U + t.

Sums over r, r' in U are split into r = r' (u terms) and r != r' (u(u - 1) terms).  Vectors are combinations of
e_t, 1_U, e_r, e_{r'}, f (f a vector on U with sum 0, |f|^2 = F2) with coefficients that may contain f_r, f_{r'}.

Scaled coordinates: X = s x, zeta_h = w_t om z_h (rho, q unscaled), as in G_z.
  Invariant part:  x = a I + b (J - I), z_h = z.   Standard part: x = al diag(f) + be f 1^T + ga 1 f^T, z = z f.
Output (pickled): the 2 x 2 invariant form Qinv(rho, z) and the standard eigenvalue lam_std, as rational functions of
(sig, s, u); see derive().
"""
import sympy as sp

u, sig, s = sp.symbols('u sigma s', positive=True)
wt, om = sp.sqrt(sig), sp.sqrt(s)
fr, frp, F2 = sp.symbols('f_r f_rp F2')


class Vec:
    """vector on U + t: dict basis -> coefficient; basis in {'t', 'one', 'er', 'erp', 'f'}"""

    def __init__(self, d=None):
        self.d = {k: sp.sympify(v) for k, v in (d or {}).items()}

    def __add__(self, o):
        d = dict(self.d)
        for k, v in o.d.items():
            d[k] = d.get(k, 0) + v
        return Vec(d)

    def __rmul__(self, c):
        return Vec({k: c * v for k, v in self.d.items()})

    def merge(self, same):
        if not same:
            return self
        d = {}
        for k, v in self.d.items():
            kk = 'er' if k == 'erp' else k
            d[kk] = d.get(kk, 0) + v.subs(frp, fr)
        return Vec(d)


def dotU(x, y, same):
    """<x_U, y_U> (Euclidean on U), x, y already merged if same."""
    G = {('one', 'one'): u, ('one', 'er'): 1, ('one', 'erp'): 1, ('one', 'f'): 0,
         ('er', 'er'): 1, ('erp', 'erp'): 1, ('er', 'erp'): 0, ('er', 'f'): fr, ('erp', 'f'): frp, ('f', 'f'): F2}
    tot = 0
    for k1, v1 in x.d.items():
        if k1 == 't':
            continue
        for k2, v2 in y.d.items():
            if k2 == 't':
                continue
            key = (k1, k2) if (k1, k2) in G else (k2, k1)
            tot += v1 * v2 * G[key]
    return tot


def sumU(x):
    """1_U . x_U"""
    S1 = {'one': u, 'er': 1, 'erp': 1, 'f': 0}
    return sum((v * S1[k] for k, v in x.d.items() if k != 't'), sp.Integer(0))


def ip(x, y, mt, mU, lam, nt, nU, same):
    """<x, M y> with M = diag(mt on t, mU on U) - lam nu nu^T, nu = nt e_t + nU 1_U."""
    x, y = x.merge(same), y.merge(same)
    xt, yt = x.d.get('t', 0), y.d.get('t', 0)
    val = mt * xt * yt + mU * dotU(x, y, same)
    nx = nt * xt + nU * sumU(x)
    ny = nt * yt + nU * sumU(y)
    return val - lam * nx * ny


def sum_rr(expr_same, expr_diff):
    """sum over r, r' in U of an expression given in the two cases (polynomials of degree <= 2 in f_r, f_r')."""
    def red(e, pairs):
        e = sp.expand(e)
        p = sp.Poly(e, fr, frp)
        tot = 0
        for (i, j), c in p.terms():
            tot += c * pairs[(i, j)]
        return tot
    same = red(expr_same, {(0, 0): u, (1, 0): 0, (2, 0): F2, (0, 1): 0, (1, 1): F2, (0, 2): F2})
    diff = red(expr_diff, {(0, 0): u * (u - 1), (1, 0): 0, (0, 1): 0, (2, 0): (u - 1) * F2, (0, 2): (u - 1) * F2,
                           (1, 1): -F2})
    return same + diff


def sum_r(expr):
    e = sp.expand(expr)
    p = sp.Poly(e, fr)
    pairs = {(0,): u, (1,): 0, (2,): F2}
    return sum((c * pairs[m] for m, c in p.terms()), sp.Integer(0))


def form(row, Xrr, zeta, rho, q, trX, wXw, wzeta, zz):
    """Phi for a structured X.  row(same): returns (psi_r, psi_rp) as Vec on U (psi_rp = row r'), with f_r, f_rp;
    Xrr(same): returns (X_{r r'}, X_{r' r}) for the case; zeta: Vec on U; trX, wXw, wzeta, zz: scalars."""
    S = u * s
    K = sig + S / 2
    # A-part: 1/2 sum_{r,r'} X_{rr'}^2 + sum_h X_hh^2
    x_same = Xrr(True)[0]
    x_diff = Xrr(False)[0]
    A = sp.Rational(1, 2) * sum_rr(x_same ** 2, x_diff ** 2) + sum_r(x_same ** 2)
    A += sp.Rational(1, 2) * (trX + rho) ** 2 + sp.Rational(3, 2) * rho ** 2 + zz
    cpl = 2 * q * (-trX + wXw / (2 * K) - (2 * K * rho + wt * wzeta) / (2 * K - sig)) + 2 * q ** 2
    # B_U: E_T = d I + (s/2) J on U
    d = s + sig + S / 2
    lamT = (s / 2) / (d * (d + u * s / 2))      # E_T^{-1} = I/d - lamT J

    def Tp(rp_is_r, psi, same):
        """T_{r'} psi on U, r' = r' (or r if rp_is_r)."""
        e = Vec({'er': 1}) if rp_is_r else Vec({'erp': 1})
        coord = dotU(psi.merge(same), e.merge(same), same)
        wpsi = om * sumU(psi)
        return sp.Rational(1, 2) * wpsi * e + (om / 2) * psi + (om * coord) * e

    def BU(same):
        psi_r, psi_rp = row(same)
        a1 = Tp(False, psi_r, same)      # T_{r'} psi_r
        a2 = Tp(True, psi_rp, same)      # T_r psi_{r'}
        return ip(a1, a2, 0, 1 / d, lamT, 0, 1, same)
    BUv = sum_rr(BU(True), BU(False))
    # cross B-term with E_0 on U + t
    d0t = sig + (sig + S) / 2
    d0U = s + (sig + S) / 2
    g0 = sig / d0t + u * s / d0U
    lam0 = 1 / (2 + g0)
    nt, nU = wt / d0t, om / d0U       # eta_0 = wt / d0; E_0^{-1} = D0^{-1} - eta0 eta0^T/(2 + g0)
    y = Vec({'t': rho}) + zeta
    wy = wt * rho + om * sumU(zeta)
    psi_r = row(True)[0]
    er = Vec({'er': 1})
    T2 = (sp.Rational(1, 2) * om * sumU(psi_r)) * Vec({'t': 1}) + (-wt / 2) * psi_r
    yr = dotU(zeta, er, True)
    T1y = (sp.Rational(1, 2) * wy) * er + (om / 2) * y + (om * yr) * er
    cross = ip(T2, T1y, 1 / d0t, 1 / d0U, lam0, nt, nU, True)
    crossv = sum_r(cross)
    return sp.expand(A + cpl - BUv - 2 * crossv)


def derive():
    a, b, q, rho, z = sp.symbols('a b q rho z')
    al, be, ga = sp.symbols('alpha beta gamma')
    # ---- invariant part: X = s (a I + b (J - I)), zeta = w_t om z 1
    def row_inv(same):
        pr = Vec({'one': s * b, 'er': s * (a - b)})
        prp = Vec({'one': s * b, 'erp': s * (a - b)})
        return pr, prp

    def Xrr_inv(same):
        return (s * a, s * a) if same else (s * b, s * b)
    zeta = Vec({'one': wt * om * z})
    Phi = form(row_inv, Xrr_inv, zeta, rho, q, trX=u * s * a, wXw=s * s * u * (a + (u - 1) * b),
               wzeta=om * wt * om * z * u, zz=u * (wt * om * z) ** 2)
    Phi = sp.expand(Phi)
    assert not Phi.has(wt) and not Phi.has(om), 'odd powers left (invariant)'
    elim = [a, b, q]
    keep = [rho, z]
    Hm = sp.hessian(Phi, elim + keep) / 2
    C = Hm[:3, :3]
    Bm = Hm[:3, 3:]
    Am = Hm[3:, 3:]
    Qinv = (Am - Bm.T * C.LUsolve(Bm)).applyfunc(sp.factor)
    # ---- standard part: X = s (al diag f + be f 1^T + ga 1 f^T), zeta = w_t om z f, rho = q = 0
    def row_std(same):
        # row r: al f_r e_r + be f_r 1 + ga f   (times s)
        pr = Vec({'er': s * al * fr, 'one': s * be * fr, 'f': s * ga})
        prp = Vec({'erp': s * al * frp, 'one': s * be * frp, 'f': s * ga})
        return pr, prp

    def Xrr_std(same):
        if same:
            return (s * (al * fr + be * fr + ga * fr),) * 2
        return (s * (be * fr + ga * frp), s * (be * frp + ga * fr))
    zeta_s = Vec({'f': wt * om * z})
    Phis = form(row_std, Xrr_std, zeta_s, 0, 0, trX=s * (al + be * 0 + ga * 0) * 0, wXw=0, wzeta=0,
                zz=(wt * om * z) ** 2 * F2)
    Phis = sp.expand(Phis)
    assert not Phis.has(wt) and not Phis.has(om), 'odd powers left (standard)'
    Phis = sp.expand(Phis / F2)
    assert not Phis.has(F2), 'not proportional to |f|^2'
    Hs = sp.hessian(Phis, [al, be, ga, z]) / 2
    Cs = Hs[:3, :3]
    Bs = Hs[:3, 3:]
    lam_std = sp.factor(Hs[3, 3] - (Bs.T * Cs.LUsolve(Bs))[0, 0])
    return dict(Qinv=Qinv, lam_std=lam_std, Cinv=C, Cstd=Cs, Phi=Phi, Phis=Phis)


if __name__ == '__main__':
    import pickle, time
    t0 = time.time()
    R = derive()
    print('derived in %.1f s' % (time.time() - t0))
    pickle.dump({k: sp.srepr(v) for k, v in R.items()}, open('equal_u.pkl', 'wb'))
    print('Qinv =', R['Qinv'])
    print('lam_std =', R['lam_std'])
    print('PASS')
