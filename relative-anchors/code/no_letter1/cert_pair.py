"""Certificate for the pair inequality (i): F_il < 0 for i != l, every j >= 1, every k >= j + 1, all weights.

Notation (homogeneous, |s| = a + b + r + e): a = x_{t_i}, b = x_{t_l} > 0, r = x_T - a - b >= 0, e = x_U >= 0,
x_T = a + b + r, c = x_T + 2e, g_i = a/(c - a), g_l = b/(c - b), G_i = a/(c - a + b), G_l = b/(c - b + a),
delta = (1 + g_i)(1 + g_l), eps = (1 + G_i)(1 + G_l), Phi(d) = sum_u x_u/(t_u + d), t_u = c - 2 x_u in [x_T, c).
  D = delta X Y - eps Z W,  X = 2 + Phi(0), Y = 2 + G_i + G_l + Phi(a + b), Z = 2 + g_i + Phi(a), W = 2 + g_l + Phi(b),
and F_il = -D/(Z W Y) (checked exactly in check_reduction.py).  So F_il < 0 iff D > 0.
(A) U empty (e = 0, Phi = 0): D is a rational function of (a, b, r); its numerator has only nonnegative coefficients
    and a positive monomial in (a, b) alone, its denominator is a product of positive-coefficient factors.
(B) U nonempty (e > 0): with psi(t) = (2/e + 1/t, (2 + G_i + G_l)/e + 1/(t + a + b), (2 + g_i)/e + 1/(t + a),
    (2 + g_l)/e + 1/(t + b)) one has (X, Y, Z, W) = sum_u x_u psi(t_u) (because sum_u x_u = e), hence
    D = sum_{u,v} x_u x_v Dk(t_u, t_v)/2 with the symmetric kernel
        Dk(t, t') = delta (X(t) Y(t') + X(t') Y(t)) - eps (Z(t) W(t') + Z(t') W(t)).
    Substituting t = x_T + 2e s1/(1 + s1), t' = x_T + 2e s2/(1 + s2) (s1, s2 >= 0 cover [x_T, c)), Dk = N/Den with
    Den a product of positive-coefficient factors, and N = N' + sum of terms (c/2) m s1^(k-1) s2^(k-1) (s1 - s2)^2
    (AM-GM repair of the negative coefficients at s1^k s2^k), where N' has only nonnegative coefficients and a
    positive monomial in (a, b, e) alone.  So Dk > 0 on [x_T, c)^2, and D > 0.
Usage: python cert_pair.py"""
import time
import flint
from flintfield import Field
from certlib import nonneg_report, positive_denominator, amgm_repair

def build_kernel():
    K = Field(['a', 'b', 'r', 'e', 's1', 's2'])
    a, b, r, e, s1, s2 = K.gens
    one = K(1)
    xT = a + b + r; c = xT + 2 * e
    gi = a / (c - a); gl = b / (c - b); Gi = a / (c - a + b); Gl = b / (c - b + a)
    dl = (1 + gi) * (1 + gl); ep = (1 + Gi) * (1 + Gl)

    def psi(t):
        return (2 * one / e + one / t, (2 + Gi + Gl) / e + one / (t + a + b), (2 + gi) / e + one / (t + a),
                (2 + gl) / e + one / (t + b))
    t1 = xT + 2 * e * s1 / (1 + s1)
    t2 = xT + 2 * e * s2 / (1 + s2)
    X1, Y1, Z1, W1 = psi(t1)
    X2, Y2, Z2, W2 = psi(t2)
    return K, dl * (X1 * Y2 + X2 * Y1) - ep * (Z1 * W2 + Z2 * W1), psi


if __name__ == '__main__':
    t0 = time.time()
    ok = True

    # (A) U empty
    KA = Field(['a', 'b', 'r'])
    a, b, r = KA.gens
    xT = a + b + r; c = xT
    gi = a / (c - a); gl = b / (c - b); Gi = a / (c - a + b); Gl = b / (c - b + a)
    dl = (1 + gi) * (1 + gl); ep = (1 + Gi) * (1 + Gl)
    D0 = dl * 2 * (2 + Gi + Gl) - ep * (2 + gi) * (2 + gl)
    nt, neg, pure = nonneg_report(D0.n, {0, 1})
    dok, dfac = positive_denominator(D0.d, {0, 1})
    good = neg == 0 and pure and dok
    ok &= good
    print(f'(A) U empty: numerator {nt} terms, negative coefficients {neg}, positive monomial in (a, b): {pure}; '
          f'denominator factors {dfac} positive: {dok} -> {"ok" if good else "FAIL"}', flush=True)

    K, Dk, _ = build_kernel()
    N = Dk.n
    rem, reps = amgm_repair(N, 4, 5)
    nt, neg, pure = nonneg_report(rem, {0, 1, 3})
    nt0, neg0, _ = nonneg_report(N, {0, 1, 3})
    dok, dfac = positive_denominator(Dk.d, {0, 1, 3})
    good = neg == 0 and pure and dok
    ok &= good
    print(f'(B) kernel Dk: numerator {nt0} terms with {neg0} negative coefficients; AM-GM repairs at {reps}; '
          f'remainder {nt} terms, negative coefficients {neg}, positive monomial in (a, b, e): {pure}; denominator '
          f'factors {dfac} positive: {dok} -> {"ok" if good else "FAIL"} [{time.time() - t0:.1f}s]', flush=True)
    print('PASS' if ok else 'FAIL')
