"""Equal U-weights: the (rho, z) block in the frame of typeI_G_u3/NOTES.md Section 4 and its limit u -> infinity.

In that frame zeta_t = c z_t w_U/|w_U| + zeta'_t with c = sqrt(1 + S/sigma).  At equal U-weights the (rho, z_t) block
is the invariant matrix Q_inv of equal_u.py after the change of variables z = c z_t/(w_t sqrt(S)), i.e.
M(u) = diag(1, c/(w_t sqrt S)) Q_inv diag(1, c/(w_t sqrt S)); its entries are rational in (sigma, S, u) (the square
roots cancel: c^2/(sigma S) = (sigma + S)/(sigma^2 S)).
1. Reproduces the floating-point values of typeI_G_u3/NOTES.md Section 4, item 3 (sigma = 0.004, S = 1:
   0.579, 0.530, 0.505, 0.489, 0.478, 0.470, 0.464, 0.460 for u = 1..8 and 0.448 for u = 12), now from the exact
   formula (smallest eigenvalue printed to 3 decimals).
2. The limit M(inf) = lim_{u -> inf} M(u) at fixed (sigma, S) exists; it is computed exactly, and M(inf) > 0 for all
   sigma, S > 0 is certified (diagonal entries and determinant are ratios of polynomials with positive coefficients).
   So at equal U-weights the smallest eigenvalue of this block has a positive limit, as conjectured there."""
import pickle
import sympy as sp
from sympy import *  # noqa: F401,F403

u, sig, s = sp.symbols('u sigma s', positive=True)
S = sp.symbols('S', positive=True)
R = {k: eval(v) for k, v in pickle.load(open('equal_u.pkl', 'rb')).items()}
Qinv = R['Qinv']
k2 = (sig + S) / (sig ** 2 * S)               # (c/(w_t sqrt S))^2
M = sp.Matrix([[Qinv[0, 0], Qinv[0, 1] * sp.sqrt(k2)], [Qinv[0, 1] * sp.sqrt(k2), Qinv[1, 1] * k2]])
M = M.subs(s, S / u)
ok = True
ref = {1: 0.579, 2: 0.530, 3: 0.505, 4: 0.489, 5: 0.478, 6: 0.470, 7: 0.464, 8: 0.460, 12: 0.448}
for uu in [1, 2, 3, 4, 5, 6, 7, 8, 12]:
    Mu = M.subs({u: uu, sig: sp.Rational(4, 1000), S: 1})
    ev = min(sp.N(e, 20) for e in Mu.eigenvals())
    good = abs(float(ev) - ref[uu]) < 1.5e-3
    ok &= good
    print(f'u = {uu}: smallest eigenvalue {float(ev):.4f} (Section 4 of typeI_G_u3 lists {ref[uu]})  {"ok" if good else "DIFFERENT"}')
# limit u -> infinity, fixed sigma, S
a = sp.limit(sp.cancel(M[0, 0]), u, sp.oo)
b2 = sp.limit(sp.cancel(Qinv[0, 1].subs(s, S / u) ** 2 * k2), u, sp.oo)      # (M[0,1])^2
c = sp.limit(sp.cancel(M[1, 1]), u, sp.oo)
a, b2, c = sp.factor(a), sp.factor(b2), sp.factor(c)
det = sp.factor(sp.together(a * c - b2))
print('M(inf)[0,0] =', a)
print('M(inf)[1,1] =', c)
print('M(inf)[0,1]^2 =', b2)
print('det M(inf) =', det)


def positive(expr):
    num, den = sp.fraction(sp.together(expr))
    for part in (num, den):
        co = sp.Poly(sp.expand(part), sig, S).coeffs()
        if not (all(x > 0 for x in co) or all(x < 0 for x in co)):
            return False
    nco = sp.Poly(sp.expand(num), sig, S).coeffs()[0]
    dco = sp.Poly(sp.expand(den), sig, S).coeffs()[0]
    return (nco > 0) == (dco > 0)


pa, pc, pd = positive(a), positive(c), positive(det)
ok &= pa and pc and pd
print(f'M(inf)[0,0] > 0: {pa};  M(inf)[1,1] > 0: {pc};  det M(inf) > 0: {pd}  (coefficient signs in (sigma, S))')
lam = min(sp.N(e) for e in sp.Matrix([[a, sp.sqrt(b2)], [sp.sqrt(b2), c]]).subs({sig: sp.Rational(4, 1000), S: 1}).eigenvals())
print(f'limit of the smallest eigenvalue at sigma = 0.004, S = 1: {float(lam):.4f}')
print('PASS' if ok else 'FAIL')
