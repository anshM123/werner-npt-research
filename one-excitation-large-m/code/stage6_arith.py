"""Stage 6: exact arithmetic of orb37 (case constants, conversion constants) and of the recursion thresholds."""
import sympy as sp

s2 = sp.sqrt(2)
out = []


def log(*a):
    s = ' '.join(str(x) for x in a); print(s); out.append(s)


def lt(a, b):
    """exact a < b for algebraic reals (sympy decides sign of the minimal polynomial root)."""
    d = sp.nsimplify(sp.simplify(b - a))
    return bool(sp.N(d, 60) > 0) and sp.simplify(d).is_positive


c163 = sp.Rational(163, 10)
# functional-bound cases
cT = (s2 + 2) * sp.sqrt(205) / 3
cS = sp.Rational(23, 2) * s2
cL = (8 * s2 + sp.Rational(16, 19)) / (sp.Rational(141, 266) * s2)
log("functional cases: small T %.6f, small S %.6f, large %.6f  (all < 16.3: %s)" %
    (float(cT), float(cS), float(cL), all(lt(c, c163) for c in (cT, cS, cL))))
log("  bracket 1-3/14-(4/19)(1+3/14) =", sp.Rational(1) - sp.Rational(3, 14) - sp.Rational(4, 19) * (1 + sp.Rational(3, 14)))
log("  rational certificates: 1711001/105000<163/10:", sp.Rational(1711001, 105000) < c163,
    " (sqrt2+2)sqrt205/3 < 1711001/105000:", lt(cT, sp.Rational(1711001, 105000)),
    " 2277/140<163/10:", sp.Rational(2277, 140) < c163, " cS<2277/140:", lt(cS, sp.Rational(2277, 140)),
    " 2288/141<163/10:", sp.Rational(2288, 141) < c163, " cL<2288/141:", lt(cL, sp.Rational(2288, 141)))
# rank-2 conversion
tau = sp.Rational(25, 4)
small = c163 * (tau + s2)
dt = s2 - (s2 + 2) * (2 / tau + s2 / tau ** 2)
num = 21 + 16 / tau + 2 / tau ** 2
large = num / dt
log("rank-2 conversion: d_tau = %s = %.6f ; small-case coeff %.6f ; large-case coeff %.6f ; both < 125: %s" %
    (sp.simplify(dt), float(dt), float(small), float(large), lt(small, 125) and lt(large, 125)))
log("  d_tau == (393 sqrt2 - 432)/625:", sp.simplify(dt - (393 * s2 - 432) / 625) == 0, " num == 14757/625:", num == sp.Rational(14757, 625),
    " 125 d_tau - num > 18/625:", lt(sp.Rational(18, 625), 125 * dt - num))
large_2pq = (21 + 16 / tau + 4 / tau ** 2) / dt
log("  robustness: with the weaker 2pq double-residual term (star21 (2) instead of orb26 (2)) large-case coeff = %.4f (< 125: %s)" %
    (float(large_2pq), lt(large_2pq, 125)))
log("rank-1 conversion coefficient (163/10) sqrt2 = %.6f < 24: %s" % (float(c163 * s2), lt(c163 * s2, 24)))


def K(B, m):
    return B * sum(s2 ** k for k in range(m - 1))


for name, B in (("B=125", sp.Integer(125)), ("B=(163/10)sqrt2", c163 * s2), ("B=24 (rounded)", sp.Integer(24)),
                ("B=max(small,large) exact", sp.Max(small, large)), ("B=4 (conjectured sharp)", sp.Integer(4))):
    ok = [m for m in range(1, 41) if lt(K(B, m), 2 ** m - 2) or sp.simplify(K(B, m) - (2 ** m - 2)) == 0]
    first = None
    for m in range(1, 41):
        if all((m2 in ok) for m2 in range(m, 41)):
            first = m; break
    log(f"{name}: K_m <= 2^m-2 holds for m in {ok[:8]}{'...' if len(ok) > 8 else ''}; holds for ALL m >= {first} (checked to 40)")
log("K16 with B=125:", sp.nsimplify(sp.expand(K(125, 16))), "= %.2f ; 125(255+127sqrt2) equal: %s ; < 111375/2: %s ; < 65534: %s"
    % (float(K(125, 16)), sp.simplify(K(125, 16) - 125 * (255 + 127 * s2)) == 0, lt(K(125, 16), sp.Rational(111375, 2)), lt(K(125, 16), 65534)))
log("K15 with B=125 = %.2f vs 2^15-2 = 32766  (fails)" % float(K(125, 15)))
log("K11 with B=24: %.2f = 744(1+sqrt2): %s ; < 1860 < 2046: %s ; K10 = %.2f vs 1022 (fails)" %
    (float(K(24, 11)), sp.simplify(K(24, 11) - 744 * (1 + s2)) == 0, lt(K(24, 11), 1860), float(K(24, 10))))
log("K10 with exact B=(163/10)sqrt2 = %.2f vs 1022 (fails); K11 = %.2f vs 2046" % (float(K(c163 * s2, 10)), float(K(c163 * s2, 11))))
log("continuation: 3*2^(m-1)+122 <= 2^(m+1)-2 iff 2^(m-1) >= 124 iff m >= 8;  sqrt2 < 3/2:", lt(s2, sp.Rational(3, 2)))
open('stage6_arith.log', 'w').write('\n'.join(out) + '\n')
