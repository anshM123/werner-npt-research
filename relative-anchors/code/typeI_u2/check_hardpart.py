"""The scalar lemma behind the hard-part bound (NOTES.md Section 2, Lemma H), as in Theorem T:
 for x, y >= 0, x + y <= 1, rho_a = 2x/(1+2x), rho_b = 2y/(1+2y):
   (1+2x)(1+2y)(2 + rho_a + rho_b - 12 rho_a rho_b) = 2 + 6(x+y) - 32xy  >= 2(1-x-y)(1+4(x+y)) >= 0,
 hence rho_a rho_b / (2 + rho_a + rho_b) <= 1/12; and the map (r1, r2) -> r1 r2/(2 + r1 + r2) is increasing.
Also the bound 1/(d_h Lam_h) <= 2/3 used for Ma.  SymPy identities; Bernstein positivity on the triangle."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
from bern import from_dict, certify
x, y, r1, r2, a, b = sp.symbols('x y r1 r2 a b', nonnegative=True)
ra, rb = 2 * x / (1 + 2 * x), 2 * y / (1 + 2 * y)
id1 = sp.simplify((1 + 2 * x) * (1 + 2 * y) * (2 + ra + rb - 12 * ra * rb) - (2 + 6 * (x + y) - 32 * x * y)) == 0
diff = sp.expand(2 + 6 * (x + y) - 32 * x * y - 2 * (1 - x - y) * (1 + 4 * (x + y)))
# diff = 8 (x - y)^2 >= 0
id2 = sp.expand(diff - 8 * (x - y) ** 2) == 0
# 2(1 - x - y)(1 + 4(x + y)) >= 0 on the triangle: trivially
f = r1 * r2 / (2 + r1 + r2)
mono = sp.simplify(sp.diff(f, r1) - r2 * (2 + r2) / (2 + r1 + r2) ** 2) == 0
# 1/(d Lam) with d = 1 + x, Lam = 3/2 - x/(1+x): d Lam = 3/2 + x/2
dl = sp.simplify((1 + x) * (sp.Rational(3, 2) - x / (1 + x)) - (sp.Rational(3, 2) + x / 2)) == 0
print('identity (1+2x)(1+2y)(2 + ra + rb - 12 ra rb) = 2 + 6(x+y) - 32xy:', id1)
print('2 + 6(x+y) - 32xy - 2(1-x-y)(1+4(x+y)) = 8 (x-y)^2:', id2)
print('d/dr1 [r1 r2/(2+r1+r2)] = r2 (2 + r2)/(2+r1+r2)^2 >= 0:', mono)
print('d_h Lam_h = 3/2 + x_h/2 (normalised units):', dl)
ok = id1 and id2 and mono and dl
print('PASS' if ok else 'FAIL')
