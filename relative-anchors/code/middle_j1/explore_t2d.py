"""Exploration (SymPy): P = det-polynomial of the T={t,t'} block in the cone variables F0, d1, d2, d3."""
import sympy as sp

al, be, F0, d1, d2, d3, p, q = sp.symbols('alpha beta F0 d1 d2 d3 p q', positive=True)
ga = al / (2 * (1 - al)); gb = be / (2 * (1 - be))
kap = (1 + ga) * (1 + gb)
G0 = 1 + (ga + gb) / 2
ap = 1 + al / (2 - 2 * al - be); bp = 1 + be / (2 - 2 * be - al)
K1 = (1 + ap) / 2; K2 = (1 + bp) / 2
F1 = F0 + d1; F2 = F0 + d2; F3 = F0 + d1 + d2 + d3
Gam = G0 + F0 / 2; Gp = K1 + F1 / 2; Gpp = K2 + F2 / 2; G3 = 1 + F3 / 2
P = 4 * Gam**2 * Gp**2 * Gpp**2 - kap**2 / 4 * Gp**2 * Gpp**2 - G3 * (ap**2 * Gam**2 * Gpp**2 + bp**2 * Gam**2 * Gp**2 - kap * ap * bp * Gam * Gp * Gpp / 2)
# triangle substitution alpha = p/(1+p+q), beta = q/(1+p+q)
sub = {al: p / (1 + p + q), be: q / (1 + p + q)}
Pp = sp.together(P.subs(sub))
num, den = sp.fraction(Pp)
num = sp.expand(num)
print('den factors:', sp.factor(den))
poly = sp.Poly(num, F0, d1, d2, d3)
neg = 0
for mon, cf in poly.terms():
    cp = sp.Poly(sp.expand(cf), p, q)
    cs = cp.coeffs()
    sgn = 'all+' if all(c > 0 for c in cs) else ('all-' if all(c < 0 for c in cs) else 'mixed')
    if sgn != 'all+':
        neg += 1
    print(mon, sgn, len(cs), 'min coef', min(cs), 'max', max(cs))
print('monomials not all-positive:', neg)
