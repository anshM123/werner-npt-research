"""Certificate for the row inequality (ii'): for every pair i != l,
    H_il := a_l F_ii + (x_T - a_i) F_il > 0,
every j >= 1, every k >= j + 1, all weights.  (With (i) this gives F_ii > sum_{l != i} |F_il|: sum H_il over l.)

Notation as in cert_pair.py, with p = a_i, q = a_l > 0, r = x_T - p - q >= 0, e = x_U >= 0, c = x_T + 2e,
alpha1 = c - p, alpha2 = c - q, beta1 = c - p + q, beta2 = c - q + p (so g_i = p/alpha1, ..., G_l = q/beta2).
Then F_ii = 2 - X (1 + g_i)^2/Z^2, F_il = -D/(Z W Y) and H_il = N_H/(Z^2 W Y) with
    N_H = 2q Z^2 W Y - q (1 + g_i)^2 X W Y - (q + r) delta X Y Z + (q + r) eps Z^2 W        (check_reduction.py).
(A) U empty (e = 0): H_il is a rational function of (p, q, r) whose numerator has only nonnegative coefficients
    and a positive monomial in (p, q) alone, and whose denominator is a product of positive-coefficient factors.
(B) U nonempty: multiply by Lambda = alpha1^2 alpha2 beta1 beta2 and put Zt = alpha1 Z, Wt = alpha2 W,
    Yt = beta1 beta2 Y.  Then Lambda N_H = 2q Zt^2 Wt Yt - q c^2 X Wt Yt - (q + r) c^2 X Yt Zt
    + (q + r)(beta1 + p)(beta2 + q) Zt^2 Wt.  Each of X, Yt, Zt, Wt is sum_u x_u f(t_u) for an explicit per-atom
    function f, and the degree-3 terms are multiplied by 1 = sum_u x_u (1/e).  Hence
        Lambda N_H = sum_{u1..u4} x_u1 x_u2 x_u3 x_u4 P(t_u1, ..., t_u4) / (24 e^4 prod_i den(t_ui)),
    den(t) = t (t + p)(t + q)(t + p + q), where P is the sum over the 24 assignments of the four letters of each
    word to the four slots of the slot polynomials e den(t) f(t).  With t_i = x_T + s_i (s_i >= 0), P is a polynomial
    in (p, q, r, e, s1, s2, s3, s4) with only nonnegative coefficients and a positive monomial in (p, q, e) alone.
    So N_H > 0, hence H_il > 0, for every finite measure of mass e on [x_T, infinity), in particular for every U.
Usage: python cert_row.py"""
import time, itertools
import flint
from flintfield import Field
from certlib import nonneg_report, positive_denominator

VARS = ('p', 'q', 'r', 'e', 's1', 's2', 's3', 's4')
WORDS_SIGN = [('2q', 'ZZWY'), ('-q c^2', 'XWYM'), ('-(q + r) c^2', 'XYZM'), ('(q + r)(beta1 + p)(beta2 + q)', 'ZZWM')]


def slot_polys(ctx, s):
    p, q, r, e = ctx.gens()[:4]
    xT = p + q + r
    al1 = q + r + 2 * e; al2 = p + r + 2 * e; be1 = 2 * q + r + 2 * e; be2 = 2 * p + r + 2 * e
    t = xT + s
    den = t * (t + p) * (t + q) * (t + p + q)
    return {'X': 2 * den + e * (t + p) * (t + q) * (t + p + q),                       # e den (2/e + 1/t)
            'Z': (2 * al1 + p) * den + e * al1 * t * (t + q) * (t + p + q),             # e den alpha1 ((2 + g_i)/e + 1/(t + p))
            'W': (2 * al2 + q) * den + e * al2 * t * (t + p) * (t + p + q),             # e den alpha2 ((2 + g_l)/e + 1/(t + q))
            'Y': (2 * be1 * be2 + p * be2 + q * be1) * den + e * be1 * be2 * t * (t + p) * (t + q),
            'M': den}                                                                    # e den (1/e)


def build_P():
    ctx = flint.fmpq_mpoly_ctx.get(VARS, 'lex')
    p, q, r, e, s1, s2, s3, s4 = ctx.gens()
    c = p + q + r + 2 * e
    be1 = 2 * q + r + 2 * e; be2 = 2 * p + r + 2 * e
    coefs = [2 * q, -q * c * c, -(q + r) * c * c, (q + r) * (be1 + p) * (be2 + q)]
    slots = [slot_polys(ctx, s) for s in (s1, s2, s3, s4)]
    P = ctx.from_dict({})
    for coef, (_, w) in zip(coefs, WORDS_SIGN):
        acc = ctx.from_dict({})
        for sg in itertools.permutations(range(4)):
            acc = acc + slots[sg[0]][w[0]] * slots[sg[1]][w[1]] * slots[sg[2]][w[2]] * slots[sg[3]][w[3]]
        P = P + coef * acc
    return ctx, P


def H_empty():
    K = Field(['p', 'q', 'r'])
    p, q, r = K.gens
    xT = p + q + r; c = xT
    gi = p / (c - p); gl = q / (c - q); Gi = p / (c - p + q); Gl = q / (c - q + p)
    dl = (1 + gi) * (1 + gl); ep = (1 + Gi) * (1 + Gl)
    X = K(2); Y = 2 + Gi + Gl; Z = 2 + gi; W = 2 + gl
    D = dl * X * Y - ep * Z * W
    Fii = 2 - X * (1 + gi) * (1 + gi) / (Z * Z)
    return q * Fii - (q + r) * D / (Z * W * Y)


if __name__ == '__main__':
    t0 = time.time()
    ok = True
    HA = H_empty()
    nt, neg, pure = nonneg_report(HA.n, {0, 1})
    dok, dfac = positive_denominator(HA.d, {0, 1})
    good = neg == 0 and pure and dok
    ok &= good
    print(f'(A) U empty: numerator {nt} terms, negative coefficients {neg}, positive monomial in (p, q): {pure}; '
          f'denominator factors {dfac} positive: {dok} -> {"ok" if good else "FAIL"}', flush=True)
    ctx, P = build_P()
    nt, neg, pure = nonneg_report(P, {0, 1, 3})
    print(f'(B) 4-kernel numerator P: {nt} terms, total degree {P.total_degree()}, degrees {P.degrees()}, '
          f'negative coefficients {neg}, positive monomial in (p, q, e): {pure} [{time.time() - t0:.1f}s]', flush=True)
    good = neg == 0 and pure
    ok &= good
    # the slot denominators e den(t) and Lambda are products of positive-coefficient polynomials
    Kd = flint.fmpq_mpoly_ctx.get(VARS, 'lex')
    p, q, r, e, s1 = Kd.gens()[:5]
    t = p + q + r + s1
    dd = e * t * (t + p) * (t + q) * (t + p + q) * (q + r + 2 * e) ** 2 * (p + r + 2 * e) * (2 * q + r + 2 * e) * (2 * p + r + 2 * e)
    dok, dfac = positive_denominator(dd, {0, 1, 3})
    ok &= dok
    print(f'    e den(t) Lambda: factors {dfac} positive: {dok}', flush=True)
    print('PASS' if ok else 'FAIL')
