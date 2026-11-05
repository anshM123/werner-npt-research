"""Closed forms of the u = 1 blocks via the structure of Phi (README, Proposition U1).
Scaled coordinates: x_t = (rho_t, zt_t) with zeta_t = w_t * om * zt_t.  Symbols: s (= s_t), sp_ (= s_t'), S (= S_T),
sg (= sigma = om^2).  Output: Gt_tt (2x2) and Gt_tt' (2x2) as rational functions; then check against u1form numerically."""
import sympy as sp, pickle
s, s2, S, sg = sp.symbols('s s2 S sigma', positive=True)
w, w2, om = sp.sqrt(s), sp.sqrt(s2), sp.sqrt(sg)
e = S + 2 * sg; K = S + sg / 2
def det(st): return (S - st / 2 + 2 * sg) * (S + st + sg / 2) - sg * st / 4
def coeffs(st, wt):
    """per-copy coefficient vectors (rho, zeta) [UNSCALED] of: alpha_t (local), A-term w_t*alpha_t, gamma0, ell"""
    d = det(st)
    a = -1 - (st * (S + st + sg) + sg * e) / (2 * d)       # zeta coef of alpha
    ap = -(om * wt * e) / d                                # rho coef of alpha
    g = wt * (S + st + sg) / d                             # zeta coef of gamma0
    gp = -om * e / d                                       # rho coef of gamma0
    c1 = -(2 * K) / (2 * K - st); c2 = -(wt * om) / (2 * K - st)
    return dict(alpha=(ap, a), A=(wt * ap, wt * a), gam=(gp, g), ell=(c1, c2), R=(1, 0))
gq = 2 * om / K
kappa = 2 - (S / 4) * gq ** 2 - S ** 2 * gq ** 2 / (8 * sg)
# global quadratic form in functionals (R, A, gam, ell):  1/2 R^2 - 1/2 A gam - S/4 gam^2 - (A + S gam)^2/(8 sg) - Lam^2/kappa
# with Lam = ell - gq/4 A - S gq/4 gam - S gq/(8 sg) (A + S gam)
fn = ['R', 'A', 'gam', 'ell']
Rv, Av, Gv, Lv = sp.symbols('Rv Av Gv Lv')
Lam = Lv - gq / 4 * Av - S * gq / 4 * Gv - S * gq / (8 * sg) * (Av + S * Gv)
glob = sp.Rational(1, 2) * Rv ** 2 - sp.Rational(1, 2) * Av * Gv - S / 4 * Gv ** 2 - (Av + S * Gv) ** 2 / (8 * sg) - Lam ** 2 / kappa
vars_ = [Rv, Av, Gv, Lv]
Q = sp.Matrix(4, 4, lambda i, j: sp.Rational(1, 2) * sp.diff(glob, vars_[i], vars_[j]))
def V(st, wt):
    c = coeffs(st, wt)
    return sp.Matrix([[c[f][0], c[f][1]] for f in fn])    # 4 x 2
Vt, Vt2 = V(s, w), V(s2, w2)
glob_off = Vt.T * Q * Vt2          # contribution to G_tt' (bilinear, x_t^T (.) x_t')
glob_diag = Vt.T * Q * Vt
# local terms at t: (3/2) rho^2 + zeta^2 - 1/4 alpha_t^2   (R^2's diagonal part is inside glob_diag)
ca = coeffs(s, w)['alpha']
loc = sp.Matrix([[sp.Rational(3, 2), 0], [0, 1]]) - sp.Rational(1, 4) * sp.Matrix([[ca[0] * ca[0], ca[0] * ca[1]], [ca[1] * ca[0], ca[1] * ca[1]]])
# pair term for ordered (t, t') and (t', t): - <M1(x_t), Ei M2(x_t')> - <M1'(x_t'), Ei' M2'(x_t)>; as bilinear in (x_t, x_t')
def pair_bilinear(st, wt, su, wu):
    """returns 2x2 B with  term = x_t^T B x_u  for the ordered pair (t, u) [coords 0, t, u]"""
    base = K - st / 2 - su / 2
    wF = [om, wt, wu]; sF = [sg, st, su]
    E = sp.Matrix(3, 3, lambda i, j: (sF[i] + base if i == j else 0) + wF[i] * wF[j] / 2)
    Ei = E.inv()
    # M1 as 3 x 2 (rows: coord 0, t, u; cols: rho_t, zeta_t)
    M1 = sp.Matrix([[0, -wu / 2], [-wu / 2, 0], [wt / 2, om / 2]])
    # M2 as 3 x 2 (rows: coord 0, t, u; cols: rho_u, zeta_u)
    M2 = sp.Matrix([[0, -wt / 2], [wu / 2, om / 2], [-wt / 2, 0]])
    return -(M1.T * Ei * M2)
Bp = pair_bilinear(s, w, s2, w2) + pair_bilinear(s2, w2, s, w).T     # both orderings, as x_t^T (.) x_t'
# total off-diagonal block (as bilinear x_t^T G_tt' x_t'): the quadratic form has 2 x_t^T G_tt' x_t', and the pair loop
# contributes both orderings once each -> G_tt' = Bp/2 ; global: x^T Q x expands to 2 x_t^T (Vt^T Q Vt2) x_t' -> G_tt' += glob_off
Goff = Bp / 2 + glob_off
Gdiag = loc + glob_diag
# scaled coordinates: zeta = w om zt
Dt = sp.diag(1, w * om); Dt2 = sp.diag(1, w2 * om)
Gt_off = (Dt * Goff * Dt2).applyfunc(lambda x: sp.factor(sp.simplify(x)))
Gt_diag = (Dt * Gdiag * Dt).applyfunc(lambda x: sp.factor(sp.simplify(x)))
pickle.dump((Gt_off, Gt_diag), open('gstruct.pkl', 'wb'))
for i in range(2):
    for j in range(2):
        print('Gt_off', i, j, ':', Gt_off[i, j], flush=True)
for i in range(2):
    for j in range(2):
        print('Gt_diag', i, j, ':', Gt_diag[i, j], flush=True)
print('closed forms written to gstruct.pkl')
