"""Negative controls for prove_conditions.py: the same pipeline must FAIL for
 (II) with tau = 7/10 (sup COL = 0.7116 > 0.7), (I) with tau = 2, (III) for Q + 1/100 (Q = 0 at the corner);
and for each a witness point with the wrong sign is exhibited (exact).  Prints PASS if all controls behave."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
from math import lcm
from fractions import Fraction as Fr
from bern import from_dict, certify
from rowbounds2 import bounds, pair_Q

s0, xi, th = sp.symbols('s0 xi theta', nonnegative=True)
a, b, c = sp.symbols('a b c', nonnegative=True)
ok = True


def cert(expr, gens, corner=None, max_boxes=3000):
    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    P = sp.Poly(sp.expand(num), *gens)
    dd = 1
    for co in P.coeffs():
        dd = lcm(dd, int(sp.denom(co)))
    dct = {m: int(co * dd) for m, co in P.terms()}
    # the denominators here are positive on the box (same factors as in prove_conditions.py)
    return certify(from_dict(dct, tuple(P.degree_list())), corner=corner, max_boxes=max_boxes)


r = bounds(2 * s0, xi, th, sp.Integer(0), one=sp.Integer(1))
res = cert(sp.Rational(7, 10) - r['COLp'], (s0, xi, th))
wit = bounds(Fr(138, 1000), Fr(1), Fr(0), Fr(0), one=Fr(1))
print('(II) tau = 7/10: certificate', 'FAILED (as it must)' if not res[0] else 'passed (BAD)', '; witness 7/10 - COLp at S = 0.138, xi = 1, theta = nu = 0:',
      float(Fr(7, 10) - wit['COLp']))
ok &= (not res[0]) and (Fr(7, 10) - wit['COLp'] < 0)
r1 = bounds(2 * s0, xi, th, sp.Integer(1), one=sp.Integer(1))
res = cert(r1['LB'] - 2 * r1['ROWm'], (s0, xi, th))
best = None
for Sn in range(1, 200, 7):
    for xn in (Fr(1, 10 ** 6), Fr(1, 2), Fr(1)):
        for tn in (Fr(0), Fr(1)):
            for ν in (Fr(0), Fr(1)):
                v = bounds(Fr(Sn, 100), xn, tn, ν, one=Fr(1))
                val = v['LB'] - 2 * max(v['ROWp'], v['ROWm'])
                if best is None or val < best[0]:
                    best = (val, Fr(Sn, 100), xn, tn, ν)
print('(I) tau = 2: certificate', 'FAILED (as it must)' if not res[0] else 'passed (BAD)', '; witness LB - 2 ROW =', float(best[0]), 'at', [str(q) for q in best[1:]])
ok &= (not res[0]) and best[0] < 0
xh = 2 * a; xk = 2 * (1 - a) * b
rh = xh / (1 + xh); rk = xk / (1 + xk)
g = rh + rk + c * (2 - rh - rk)
Q = pair_Q(xh, xk, g, one=sp.Integer(1))
res = cert(-(Q + sp.Rational(1, 100)), (a, b, c), corner=(0, 0, 0))
w0 = -(pair_Q(Fr(0), Fr(0), Fr(0), one=Fr(1)) + Fr(1, 100))
print('(III) Q + 1/100: certificate', 'FAILED (as it must)' if not res[0] else 'passed (BAD)', '; witness value at the corner:', float(w0))
ok &= (not res[0]) and w0 < 0
print('PASS' if ok else 'FAIL')
