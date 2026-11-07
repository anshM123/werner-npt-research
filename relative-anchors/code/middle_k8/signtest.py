"""The sign test of the certificates: is a rational function f of q (ufield.RF) > 0 for every q > 0?

positive_on_halfline(f) returns (ok, how, witness):
 1. numerator and denominator have only nonnegative coefficients (and are not zero): ok;
 2. otherwise, a rational q > 0 at which f is <= 0 or undefined, found by exact evaluation at candidate points
    (a fixed grid first, then points next to floating-point roots, which only propose candidates): not ok, with the
    witness;
 3. otherwise the exact test of ufield.rf_positive (Sturm counts of the roots in (0, inf)).
Step 2 only shortens the negative cases (a Sturm count of a polynomial of degree about 150 takes minutes); it never
turns a negative answer into a positive one."""
from fractions import Fraction as Fr
from ufield import rf_positive

GRID = [Fr(1, 10 ** 6), Fr(1, 10 ** 4), Fr(1, 1000), Fr(1, 100), Fr(1, 30), Fr(1, 10), Fr(1, 3), Fr(1), Fr(3),
        Fr(10), Fr(30), Fr(100), Fr(1000), Fr(10 ** 4), Fr(10 ** 6), Fr(1, 10 ** 8), Fr(10 ** 8)]


def _bad(f, x):
    try:
        return f(x) <= 0
    except ZeroDivisionError:
        return True


def find_negative_point(f):
    """a rational q > 0 at which f is <= 0 or undefined, or None"""
    for x in GRID:
        if _bad(f, x):
            return x
    roots = []
    try:
        import numpy as np
        for p in (f.n, f.d):
            cs = [Fr(int(c.p), int(c.q)) for c in p.coeffs()]
            if len(cs) > 1:
                big = max(abs(c) for c in cs)
                cs = [float(c / big) for c in cs]
                for r in np.roots(cs[::-1]):
                    if r.real > 0 and abs(r.imag) <= 1e-7 * abs(r):
                        roots.append(float(r.real))
    except Exception:
        roots = []
    roots = sorted(roots)
    cands = []
    for i, r in enumerate(roots):
        cands += [Fr(r * 0.999).limit_denominator(10 ** 9), Fr(r * 1.001).limit_denominator(10 ** 9)]
        if i + 1 < len(roots):
            cands.append(Fr((r + roots[i + 1]) / 2).limit_denominator(10 ** 12))
    for x in cands:
        if x > 0 and _bad(f, x):
            return x
    return None


def positive_on_halfline(f):
    nn = (not f.n.is_zero()) and all(c >= 0 for c in f.n.coeffs())
    dd = (not f.d.is_zero()) and all(c >= 0 for c in f.d.coeffs())
    if nn and dd:
        return True, (f'num ({f.n.degree()}, nonnegative coefficients); den ({f.d.degree()}, nonnegative '
                      f'coefficients)'), None
    x = find_negative_point(f)
    if x is not None:
        return False, f'nonpositive at q = {x} (exact evaluation)', x
    ok, how = rf_positive(f)
    return ok, how, None
