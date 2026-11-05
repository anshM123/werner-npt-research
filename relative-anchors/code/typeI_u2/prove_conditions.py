"""Exact Bernstein certificates (integer arithmetic, bern.py) for the conditions built by build_conditions.py.
For each condition num/den: every denominator factor is shown to have a constant strict sign on [0,1]^3 (all
Bernstein coefficients of one strict sign), and then sign(den) * num > 0 on [0,1]^3 (all Bernstein coefficients
positive after dyadic subdivision; for (III) the corner rule at the single vertex (0,0,0), where -Q = 0).
Prints PASS."""
import sys, os, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
from math import lcm
from bern import from_dict, certify, bernstein_on_box

t00 = time.time()
data = pickle.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'conditions.pkl'), 'rb'))
ok_all = True


def int_dict(P):
    den = 1
    for co in P.coeffs():
        den = lcm(den, int(sp.denom(co)))
    return {m: int(co * den) for m, co in P.terms()}, tuple(P.degree_list())


def const_sign(P):
    """+1 / -1 if all Bernstein coefficients on the root box (or after subdivision) have that strict sign."""
    dct, degs = int_dict(P)
    C = from_dict(dct, degs)
    okp, _ = certify(C, max_boxes=20000)
    if okp:
        return 1
    Cn = from_dict({m: -v for m, v in dct.items()}, degs)
    okn, _ = certify(Cn, max_boxes=20000)
    return -1 if okn else 0


from fractions import Fraction as Fr
import random
from rowbounds2 import bounds, pair_Q
random.seed(3)


def direct(name, pt):
    """the condition evaluated directly with Fractions (no SymPy) at a rational point of [0,1]^3"""
    if name == 'III':
        a_, b_, c_ = pt
        xh = 2 * a_; xk = 2 * (1 - a_) * b_
        rh = xh / (1 + xh); rk = xk / (1 + xk)
        return -pair_Q(xh, xk, rh + rk + c_ * (2 - rh - rk), one=Fr(1))
    kind, sgn_, nu = name.split('_')
    r = bounds(2 * pt[0], pt[1], pt[2], Fr(int(nu[2:])), one=Fr(1))
    if kind == 'I':
        return r['LB'] - Fr(4, 5) * r['ROW' + sgn_]
    if kind == 'Iq':
        return r['LB'] - Fr(4, 5) * r['ROW' + sgn_] - Fr(1, 10)
    if kind == 'IIq':
        return Fr(3, 4) - r['COL' + sgn_]
    return Fr(4, 5) - r['COL' + sgn_]


for name, (num, den) in data.items():
    t0 = time.time()
    gens = num.gens
    # consistency of the SymPy rational function with the direct Fraction evaluation
    for rep in range(4):
        pt = tuple(Fr(random.randint(1, 99), 100) for _ in range(3))
        sub = dict(zip(gens, [sp.Rational(q.numerator, q.denominator) for q in pt]))
        val = sp.Rational(num.as_expr().subs(sub)) / (den[0] * sp.prod([f.subs(sub) ** e for f, e in den[1]]))
        if Fr(int(sp.numer(val)), int(sp.denom(val))) != direct(name, pt):
            ok_all = False
            print(f"  {name}: SymPy form differs from direct evaluation at {pt}", flush=True)
    sgn = 1 if den[0] > 0 else -1
    fac_ok = True
    for f, e in den[1]:
        Pf = sp.Poly(f, *gens)
        sf = const_sign(Pf)
        if sf == 0:
            fac_ok = False
            print(f"  {name}: denominator factor of degree {Pf.total_degree()} has no certified constant sign", flush=True)
        if e % 2 == 1:
            sgn *= sf if sf != 0 else 1
    dct, degs = int_dict(num)
    if sgn < 0:
        dct = {m: -v for m, v in dct.items()}
    C = from_dict(dct, degs)
    corner = (0, 0, 0) if name == 'III' else None
    ok, st = certify(C, corner=corner, max_boxes=400000)
    good = ok and fac_ok
    ok_all &= good
    print(f"{name}: denominator sign {sgn:+d} (factors certified: {fac_ok}); numerator {'certified' if ok else 'FAILED'} {st} "
          f"({time.time() - t0:.1f} s)", flush=True)
print(f"total {time.time() - t00:.1f} s")
print("PASS" if ok_all else "FAIL")
