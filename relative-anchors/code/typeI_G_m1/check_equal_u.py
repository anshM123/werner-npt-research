"""Equal U-weights, every u: exact validation of the symbolic reduction (equal_u.py) against the exact engine, and the
sign certificates that prove Conjecture M on this family for every u.

Notation: sigma = s_t, s = common U-weight squared, S = u s.  From equal_u.py:
  Qinv = the 2 x 2 matrix of G_z restricted to the S_u-invariant letter-2 data (rho, z_h = z for all h), i.e.
         Qinv[0,0] = G_z[rho,rho], Qinv[0,1] = sum_h G_z[rho,h], Qinv[1,1] = sum_{h,h'} G_z[h,h'];
  lam_std = G_z[h,h] - G_z[h,h'] (h != h'), the eigenvalue of G_z on the standard representation (u >= 3).
Part 1 (validation): exact equality with G_z from gzexact.py (exact elimination of the whole block over Q) for
u = 1, ..., 12 at three rational weight pairs each (lam_std for u >= 3).
Part 2 (certificates): with u = u0 + x, every factor of the numerator and denominator of the following rational
functions is a polynomial in (sigma, s, x) with only positive coefficients (up to an explicit overall sign):
  (Z1) -Qinv[0,1] > 0                        (G_z[rho,h] < 0), u >= 1;
  (D1) Qinv[0,0] v0 + Qinv[0,1] v1 > 0         ((G_z v)_rho > 0), u >= 1;
  (D2) Qinv[0,1] v0 + Qinv[1,1] v1 > 0         ((G_z v)_h > 0, times u), u >= 1;
  (Z2) lam_std - Qinv[1,1]/u > 0              (G_z[h,h'] < 0, times -u), u >= 3;
with v0 = 4 sigma (sigma + S), v1 = 3 (2 sigma + S).  Since x >= 0, sigma, s > 0 these hold for every real u >= u0, in
particular every integer u.  (For u = 2 the off-diagonal z-entry is covered by check_zpattern.py.)"""
import sys, os, pickle, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
import flint
from sympy import *  # noqa: F401,F403  (namespace for srepr)
from gzexact import Gz_exact

Q = flint.fmpq
u, sig, s = sp.symbols('u sigma s', positive=True)
x = sp.symbols('x', nonnegative=True)
R = {k: eval(v) for k, v in pickle.load(open('equal_u.pkl', 'rb')).items()}
Qinv, lam = R['Qinv'], R['lam_std']
ok = True

# --- integrality of exponents (the odd powers of w_t, omega cancelled) in the two forms
for key in ('Phi', 'Phis'):
    for p in R[key].atoms(sp.Pow):
        if p.base in (sig, s) and not p.exp.is_integer:
            print(f'{key}: non-integer exponent {p}')
            ok = False
print('forms are polynomial in (sigma, s) (no square roots):', ok)

# --- Part 1: validation against the exact engine
def ev(expr, st, ss, uu):
    return sp.cancel(expr.subs({sig: st, s: ss, u: uu}))


pts = [(Q(3, 4), Q(1, 2)), (Q(2), Q(1, 3)), (Q(1, 7), Q(5, 2))]
nval = 0
t0 = time.time()
for uu in range(1, 13):
    for wt, om in pts:
        G = Gz_exact(wt, [om] * uu)
        st, ss = sp.Rational(int((wt * wt).p), int((wt * wt).q)), sp.Rational(int((om * om).p), int((om * om).q))
        n = uu + 1
        val = lambda q: sp.Rational(int(q.p), int(q.q))
        g00 = val(G[0, 0])
        g0h = sum(val(G[0, j]) for j in range(1, n))
        ghh = sum(val(G[i, j]) for i in range(1, n) for j in range(1, n))
        e = [ev(Qinv[0, 0], st, ss, uu) - g00, ev(Qinv[0, 1], st, ss, uu) - g0h, ev(Qinv[1, 1], st, ss, uu) - ghh]
        if uu >= 3:
            e.append(ev(lam, st, ss, uu) - (val(G[1, 1]) - val(G[1, 2])))
        good = all(t == 0 for t in e) and all(isinstance(t, sp.Rational) or t == 0 for t in e)
        ok &= good
        nval += 1
        if not good:
            print(f'   MISMATCH u = {uu}, w_t = {wt}, om = {om}: {e}')
print(f'Part 1: symbolic invariant/standard reduction = exact elimination at {nval} points (u = 1..12): '
      f'{"all equal" if ok else "MISMATCH"} ({time.time() - t0:.1f} s)')


# --- Part 2: sign certificates
def sign_cert(expr, u0, name):
    """expr rational in (sigma, s, u); returns True if expr > 0 for sigma, s > 0, u >= u0 by positive coefficients."""
    num, den = sp.fraction(sp.together(expr))
    sgn = 1
    report = []
    for part, which in ((num, 'num'), (den, 'den')):
        c, facs = sp.factor_list(sp.expand(part))
        if c < 0:
            sgn = -sgn
        for f, e in facs:
            p = sp.Poly(sp.expand(f.subs(u, u0 + x)), sig, s, x)
            co = p.coeffs()
            if all(cc > 0 for cc in co):
                report.append((which, p.total_degree(), len(co), '+', e))
            elif all(cc < 0 for cc in co):
                report.append((which, p.total_degree(), len(co), '-', e))
                if e % 2 == 1:
                    sgn = -sgn
            else:
                print(f'   {name}: factor with mixed coefficients: {f}')
                return False, report
    good = sgn > 0
    print(f'{name}: u >= {u0}: every factor sign-definite; overall sign {"+" if good else "-"}; '
          f'factors (part, degree, terms, sign, exponent): {report}')
    return good, report


S = u * s
v0 = 4 * sig * (sig + S)
v1 = 3 * (2 * sig + S)
t0 = time.time()
r1, _ = sign_cert(-Qinv[0, 1], 1, '(Z1) -G_z[rho,h]')
r2, _ = sign_cert(Qinv[0, 0] * v0 + Qinv[0, 1] * v1, 1, '(D1) (G_z v)_rho')
r3, _ = sign_cert(Qinv[0, 1] * v0 + Qinv[1, 1] * v1, 1, '(D2) u (G_z v)_h')
r4, _ = sign_cert(lam - Qinv[1, 1] / u, 3, '(Z2) -u G_z[h,hp]')
r5, _ = sign_cert(lam, 3, '(S) lam_std')
ok &= r1 and r2 and r3 and r4 and r5
print(f'Part 2 done ({time.time() - t0:.1f} s)')
print('PASS' if ok else 'FAIL')
