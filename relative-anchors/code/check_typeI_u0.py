"""EXACT check of Proposition U0 (relative-anchors/README.md): the type (I) block of (R_inf)_{k-1} in which EVERY copy carries the letter 2
once (T = all copies, no letter 1 anywhere; states R_t, t < k, and Q) is positive definite for every k >= 2.
Proof skeleton checked here:
 (1) closed form (in units of 2Y): 2 on the R_t diagonal, 1/2 + kappa_tt'/4 off the diagonal, -c_t = -2S/(2S - s_t)
     between R_t and Q, 2 at Q  [exact, against the builder, k = 2..5];
 (2) the Schur complement of Q equals F with F_tt = D_t = 2(S - s_t)(3S - s_t)/(2S - s_t)^2 and
     F_tt' = beta(s_t, s_t') = kappa/4 - (d_t + d_t')/2 - d_t d_t'/2, d_t = s_t/(2S - s_t)  [symbolic identity];
 (3) with c = S - a - b >= 0:  beta(a,b) = -a b P1/Den,  P1, Den with nonnegative coefficients  [symbolic];
 (4) 2 b (3S - a)/(2S - a)^2 + beta(a,b) = b P2/Den2 with P2, Den2 with nonnegative coefficients and P2 > 0  [symbolic];
 (5) hence sum_{t' != t} |F_tt'| < (S - s_t) 2(3S - s_t)/(2S - s_t)^2 = D_t: F is strictly diagonally dominant with
     positive diagonal, so F > 0 and the block is positive definite.  Also checked at rational points (k = 2..8).
Prints PASS/FAIL."""
import random
import sympy as sp
import flint
from ylib import Builder, point_ring, rinf_states, wstate, Nv

ok = True
a, b, c = sp.symbols('a b c', nonnegative=True)
S = a + b + c


def kap_(x, y, SS):
    return (SS * (x + y) + x**2 + y**2) / (SS**2 + SS * (x + y) / 2 + x * y - (x**2 + y**2) / 2)


# (2) Schur complement identity, symbolic (entries only depend on the pair and on S)
x, y, SS = sp.symbols('x y SS', positive=True)
cx, cy = 2 * SS / (2 * SS - x), 2 * SS / (2 * SS - y)
dx, dy = x / (2 * SS - x), y / (2 * SS - y)
diag_ok = sp.simplify((2 - cx * cx / 2) - 2 * (SS - x) * (3 * SS - x) / (2 * SS - x)**2) == 0
off_ok = sp.simplify((sp.Rational(1, 2) + kap_(x, y, SS) / 4 - cx * cy / 2) - (kap_(x, y, SS) / 4 - (dx + dy) / 2 - dx * dy / 2)) == 0
print('(2) Schur complement of Q: diagonal formula', diag_ok, '; off-diagonal formula', off_ok)
ok &= diag_ok and off_ok
# (3), (4)
beta = sp.together(kap_(a, b, S) / 4 - (a / (2 * S - a) + b / (2 * S - b)) / 2 - (a / (2 * S - a)) * (b / (2 * S - b)) / 2)
P1 = 4 * a**2 + 7 * a * b + 7 * a * c + 4 * b**2 + 7 * b * c + 3 * c**2
Den = (a + 2 * b + 2 * c) * (2 * a + b + 2 * c) * (2 * a**2 + 8 * a * b + 5 * a * c + 2 * b**2 + 5 * b * c + 2 * c**2)
id3 = sp.simplify(beta + a * b * P1 / Den) == 0
pos3 = all(v >= 0 for v in sp.Poly(sp.expand(P1), a, b, c).coeffs()) and all(v >= 0 for v in sp.Poly(sp.expand(Den), a, b, c).coeffs())
print('(3) beta = -a b P1/Den:', id3, '; P1, Den nonnegative coefficients:', pos3)
expr = sp.together(2 * b * (3 * S - a) / (2 * S - a)**2 + beta)
n2, d2 = sp.fraction(sp.factor(expr))
P2n = sp.Poly(sp.expand(n2), a, b, c); P2d = sp.Poly(sp.expand(d2), a, b, c)
sgn = 1 if P2d.LC() > 0 else -1
pos4 = all(sgn * v >= 0 for v in P2n.coeffs()) and all(sgn * v >= 0 for v in P2d.coeffs())
q4, r4 = sp.div(sp.expand(n2), b, a, b, c)
has_const = sp.Poly(sp.expand(q4), a, b, c).coeffs()
print('(4) bound numerator and denominator nonnegative coefficients:', pos4, '; numerator = b * (poly with',
      'positive pure powers a^4, b^4, c^4):', r4 == 0 and all(sgn * sp.Poly(q4, a, b, c).coeff_monomial(m) > 0 for m in (a**4, b**4, c**4)))
ok &= id3 and pos3 and pos4 and r4 == 0
# (1) closed form against the builder, and (5) at rational points
Q_ = flint.fmpq
random.seed(23)
for k in range(2, 9):
    w = [Q_(random.randint(1, 15), random.randint(1, 15)) for _ in range(k)]
    if k % 2 == 0:
        w[0] = Q_(50)                      # one dominant weight
    s = [t * t for t in w]
    Ssum = sum(s, Q_(0))
    if k <= 5:
        Bd = Builder(point_ring(s))
        key = (0,) + (1,) * k
        sts = [X for X in rinf_states(k, k - 1) if wstate(*X) == key]

        def mu(X):
            v = Q_(1)
            for t, n in enumerate(Nv(*X)):
                v *= w[t] ** n
            return v
        lab = ['Q' if sum(A_) == 0 else A_.index(2) for (A_, y_) in sts]
        Y2 = {(lab[i], lab[j]): 2 * Bd.Y(sts[i], sts[j]) / (mu(sts[i]) * mu(sts[j])) for i in range(len(sts)) for j in range(len(sts))}
        cf = Y2[('Q', 'Q')] == 2
        for t in range(k):
            cf &= Y2[(t, t)] == 2 and Y2[(t, 'Q')] == -2 * Ssum / (2 * Ssum - s[t])
            for u in range(k):
                if u != t:
                    cf &= Y2[(t, u)] == Q_(1, 2) + kap_(s[t], s[u], Ssum) / 4
        ok &= cf
    else:
        cf = 'n/a (k > 5)'
    # F and strict diagonal dominance, exact
    dd = True
    for t in range(k):
        Dt = 2 * (Ssum - s[t]) * (3 * Ssum - s[t]) / (2 * Ssum - s[t])**2
        off = Q_(0)
        for u in range(k):
            if u != t:
                bt = kap_(s[t], s[u], Ssum) / 4 - (s[t] / (2 * Ssum - s[t]) + s[u] / (2 * Ssum - s[u])) / 2 \
                    - (s[t] / (2 * Ssum - s[t])) * (s[u] / (2 * Ssum - s[u])) / 2
                dd &= bt < 0
                off += -bt
        dd &= off < Dt
    ok &= dd
    print(f'k={k}: closed form vs builder {cf}; F strictly diagonally dominant with negative off-diagonal: {dd}')
print('PASS' if ok else 'FAIL')
