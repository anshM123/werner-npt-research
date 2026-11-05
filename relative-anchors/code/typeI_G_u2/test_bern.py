"""Unit tests of bern.py (exact): vertex coefficients = values, interior reconstruction, a positive polynomial that
needs subdivision, a polynomial with a negative value (must fail), the corner rule (and failure without it)."""
from fractions import Fraction as Fr
from math import comb
import random
from bern import from_dict, bernstein_on_box, certify

random.seed(1)
degs = (3, 4)
dct = {(i, j): random.randint(-9, 9) for i in range(4) for j in range(5)}
C = from_dict(dct, degs)
def ev(x, y): return sum(Fr(c) * x**i * y**j for (i, j), c in dct.items())
box = ((2, 1), (3, 5))
B = bernstein_on_box(C, box)
x0, x1 = Fr(1, 4), Fr(2, 4); y0, y1 = Fr(5, 8), Fr(6, 8)
sc = Fr(2) ** (2 * 3 + 3 * 4)
ok1 = B[0, 0] == ev(x0, y0) * sc and B[3, 4] == ev(x1, y1) * sc and B[3, 0] == ev(x1, y0) * sc and B[0, 4] == ev(x0, y1) * sc
X, Y = Fr(1, 3), Fr(2, 7)
tx, ty = (X - x0) / (x1 - x0), (Y - y0) / (y1 - y0)
val = sum(Fr(B[i, j], comb(3, i) * comb(4, j)) * comb(3, i) * tx**i * (1 - tx)**(3 - i) * comb(4, j) * ty**j * (1 - ty)**(4 - j)
          for i in range(4) for j in range(5)) / sc
ok1 &= val == ev(X, Y)
print('vertex values and interior reconstruction on a sub-box:', ok1)
r1 = certify(from_dict({(0,): 1009, (1,): -6000, (2,): 9000}, (2,)))      # 9000 (x - 1/3)^2 + 9 > 0
print('9000 (x - 1/3)^2 + 9 certified:', r1)
r2 = certify(from_dict({(0,): 991, (1,): -6000, (2,): 9000}, (2,)), max_boxes=2000)   # negative near 1/3
print('9000 (x - 1/3)^2 - 9 rejected:', r2[0] is False, r2[1].get('reason'))
C3 = from_dict({(1, 0): 1, (0, 1): 1, (1, 1): 1, (2, 1): -1}, (2, 1))   # x + y + xy - x^2 y, zero only at (0,0)
r3 = certify(C3, corner=(0, 0))
r4 = certify(C3, max_boxes=500)
print('corner rule:', r3, '; without it rejected:', r4[0] is False)
ok = ok1 and r1[0] and (r2[0] is False) and r3[0] and (r4[0] is False)
print('PASS' if ok else 'FAIL')
