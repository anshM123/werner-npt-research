"""Sign certificates for the weighted block Gershgorin conditions at equal U-weights (u >= 3, every m >= 2).

Input: the closed forms of equal_m.py (validated exactly by check_equal_m.py).  Orthant variables:
  (A), off-diagonal blocks:  s_t = a, s_t' = b, sigma = a + b + r (r = weight of the other T-copies, r >= 0), U-weight s,
                             u = 3 + x;
  (B), diagonal blocks:      s_t = a, sigma = a + r (r > 0 since m >= 2), U-weight s, u = 3 + x.
Frame: invariant part (rho_t, zhat_t), zhat_t = kappa_t z_t, kappa_t^2 = s_t S sigma/(sigma + S), S = u s; standard part
zhat = w_t om |f| z.  With H the Gershgorin constant (default 1/2) the conditions are
  (A_inv)  H^2 s_t s_t'/sigma^2 - ||Ghat_tt'||_F^2 > 0      (A_std)  H^2 s_t s_t'/sigma^2 - Ghat_tt'^2 > 0
  (B_inv)  N = Ghat_tt - H (1 - theta) I > 0, i.e. N_00 > 0 and kappa^2 det N > 0      (B_std)  Ghat_tt - H(1-theta) > 0
(theta = s_t/sigma, 1 - theta = r/sigma).  Each is a rational function of (a, b, r, s, x); the certificate: the reduced
numerator has only positive coefficients (with a monomial free of x and of r), and every irreducible factor of the
reduced denominator has only positive coefficients.  Then the condition holds for all a, b, s > 0, r > 0 (r >= 0 for
(A)), x >= 0, i.e. for every real u >= 3.
Usage: python prove_em.py [H_num H_den]"""
import sys, os, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from emlib import load
from flintfield import Field, F as FE

Q = flint.fmpq
R = load()
VN = ('a', 'b', 'r', 's', 'x')
Fd = Field(list(VN))
ctx = Fd.ctx
a, b, r, s, x = ctx.gens()
one = ctx.from_dict({(0,) * 5: 1})


def comp(nd, sig, s1, s2):
    n, d = nd
    args = (3 * one + x, sig, s, s1, s2)
    return FE(Fd, n.compose(*args, ctx=ctx), d.compose(*args, ctx=ctx))


def lift(p):
    return FE(Fd, p, one, reduce=False)


def cert(name, val, need_free_of=('x',)):
    """val: rational function; returns True if certified positive on the domain."""
    t0 = time.time()
    num, den = val.n, val.d
    # sign normalisation: den has positive leading coefficient by construction of the field
    ok = True
    cn = list(num.coeffs())
    neg = sum(1 for c in cn if c < 0)
    pos = sum(1 for c in cn if c > 0)
    sgn_num = 1 if neg == 0 else (-1 if pos == 0 else 0)
    # a monomial free of the variables in need_free_of (so the polynomial is positive when they vanish)
    idx = [VN.index(v) for v in need_free_of]
    free = any(all(m[i] == 0 for i in idx) for m in num.monoms())
    c, fl = den.factor()
    den_ok = True
    nfac = len(fl)
    sgn_den = 1 if c > 0 else -1
    for f, e in fl:
        fc = list(f.coeffs())
        # each factor: one-signed coefficients and a monomial free of need_free_of (nonzero on the closed domain)
        if not any(all(m[i] == 0 for i in idx) for m in f.monoms()):
            den_ok = False
        if all(q > 0 for q in fc):
            continue
        if all(q < 0 for q in fc):
            if e % 2:
                sgn_den = -sgn_den
            continue
        den_ok = False
    ok = sgn_num == 1 and free and den_ok and sgn_den == 1
    degs = [max((m[i] for m in num.monoms()), default=0) for i in range(5)]
    print('%-10s numerator: %7d terms, degrees %s, negative coefficients %d, positive %d; monomial free of %s: %s; '
          'denominator factors one-signed: %s (sign %+d) -> %s [%.1fs]'
          % (name, len(num), degs, neg, pos, need_free_of, free, den_ok, sgn_den, 'CERTIFIED' if ok else 'not certified',
             time.time() - t0), flush=True)
    return ok, num


def conditions(H):
    """the five rational functions of (a, b, r, s, x) whose positivity is the weighted block Gershgorin criterion"""
    Ginv = R['inv']
    HQ = Fd(int(H.p)) / Fd(int(H.q))
    # (A): off-diagonal blocks
    sigA = a + b + r
    g = [[comp(Ginv[i][2 + j], sigA, a, b) for j in range(2)] for i in range(2)]
    gs = comp(R['std'][0][1], sigA, a, b)
    SA = lift((3 * one + x) * s)
    sgA = lift(sigA)
    kt2 = lift(a) * SA * sgA / (sgA + SA)
    kp2 = lift(b) * SA * sgA / (sgA + SA)
    target = HQ * HQ * lift(a) * lift(b) / (sgA * sgA)
    FA = target - (g[0][0] * g[0][0] + g[0][1] * g[0][1] / kp2 + g[1][0] * g[1][0] / kt2
                   + g[1][1] * g[1][1] / (kt2 * kp2))
    FAs = target - gs * gs / (lift(a) * lift(b) * lift(s) * lift(s))
    # (B): diagonal blocks
    sigB = a + r
    gd = [[comp(Ginv[i][j], sigB, a, 0 * one) for j in range(2)] for i in range(2)]
    gsd = comp(R['std'][0][0], sigB, a, 0 * one)
    SB = lift((3 * one + x) * s)
    sgB = lift(sigB)
    k2 = lift(a) * SB * sgB / (sgB + SB)
    om_th = HQ * lift(r) / sgB           # H (1 - theta)
    N00 = gd[0][0] - om_th
    detk = (gd[0][0] - om_th) * (gd[1][1] - om_th * k2) - gd[0][1] * gd[1][0]
    NBs = gsd / (lift(a) * lift(s)) - om_th
    return dict(A_inv=FA, A_std=FAs, B_inv_00=N00, B_inv_det=detk, B_std=NBs)


FREE = dict(A_inv=('x', 'r'), A_std=('x', 'r'), B_inv_00=('x',), B_inv_det=('x',), B_std=('x',))


def run(H, verbose=True):
    t0 = time.time()
    conds = conditions(H)
    res = {}
    for name in ('A_inv', 'A_std', 'B_inv_00', 'B_inv_det', 'B_std'):
        res[name] = cert(name, conds[name], FREE[name])[0]
    allok = all(res.values())
    print('H = %s: all five conditions certified: %s [%.1fs]' % (H, allok, time.time() - t0), flush=True)
    return allok, res, conds


if __name__ == '__main__':
    H = Q(int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else Q(1, 2)
    allok, res, conds = run(H)
    print('PASS' if allok else 'FAIL')
