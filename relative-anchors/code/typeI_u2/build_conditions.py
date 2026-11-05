"""Build the per-coordinate conditions (I), (II) of NOTES.md Section 4 and the pair condition (III) as rational
functions (SymPy), split into numerator / denominator, and save them (conditions.pkl).
 (I)   LB - (4/5) ROW_s > 0          (s = p, m: the two signs of Mca_hh), nu in {0, 1}
 (II)  4/5 - COL_s > 0               nu in {0, 1}
 (I_q) LB - (4/5) ROW_s - 1/10 > 0, (II_q) 3/4 - COL_s > 0   (quantitative versions, Corollary Q)
 Variables: S = 2 s0 (s0 in [0, 1]), xi in [0, 1], theta in [0, 1].
 (III) -Q(x_h, x_k, g) >= 0  with x_h = 2 a, x_k = 2 (1 - a) b, g = rho_h + rho_k + c (2 - rho_h - rho_k), (a, b, c) in [0,1]^3.
"""
import sys, os, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
from rowbounds2 import bounds, pair_Q

s0, xi, th = sp.symbols('s0 xi theta', nonnegative=True)
a, b, c = sp.symbols('a b c', nonnegative=True)
TAU = sp.Rational(4, 5)
t0 = time.time()
out = {}
for nu in (0, 1):
    r = bounds(2 * s0, xi, th, sp.Integer(nu), one=sp.Integer(1))
    conds = {f'I_p_nu{nu}': r['LB'] - TAU * r['ROWp'], f'I_m_nu{nu}': r['LB'] - TAU * r['ROWm'],
             f'II_p_nu{nu}': TAU - r['COLp'], f'II_m_nu{nu}': TAU - r['COLm'],
             # quantitative versions (Corollary Q): margins 1/10 and 1/20
             f'Iq_p_nu{nu}': r['LB'] - TAU * r['ROWp'] - sp.Rational(1, 10), f'Iq_m_nu{nu}': r['LB'] - TAU * r['ROWm'] - sp.Rational(1, 10),
             f'IIq_p_nu{nu}': sp.Rational(3, 4) - r['COLp'], f'IIq_m_nu{nu}': sp.Rational(3, 4) - r['COLm']}
    for name, ex in conds.items():
        num, den = sp.fraction(sp.cancel(sp.together(ex)))
        num = sp.Poly(sp.expand(num), s0, xi, th)
        den = sp.factor_list(sp.expand(den))
        out[name] = (num, den)
        print(f"{name}: numerator degrees {num.degree_list()}, {len(num.terms())} terms; denominator factors "
              f"{[(sp.Poly(f, s0, xi, th).total_degree(), e) for f, e in den[1]]} const {den[0]} ({time.time() - t0:.0f} s)", flush=True)
xh = 2 * a
xk = 2 * (1 - a) * b
rh = xh / (1 + xh); rk = xk / (1 + xk)
g = rh + rk + c * (2 - rh - rk)
Q = pair_Q(xh, xk, g, one=sp.Integer(1))
num, den = sp.fraction(sp.cancel(sp.together(-Q)))
num = sp.Poly(sp.expand(num), a, b, c)
den = sp.factor_list(sp.expand(den))
out['III'] = (num, den)
print(f"III: numerator degrees {num.degree_list()}, {len(num.terms())} terms; denominator factors "
      f"{[(sp.Poly(f, a, b, c).total_degree(), e) for f, e in den[1]]} const {den[0]} ({time.time() - t0:.0f} s)", flush=True)
pickle.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'conditions.pkl'), 'wb'))
print("saved conditions.pkl")
print("PASS")
