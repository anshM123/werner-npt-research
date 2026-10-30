"""Check the integration-by-parts recursion behind the all-n polynomial lemma.
M_r(x) = int_0^1 y^k (1-y)^(m-1) (x-y)^r dy.
General identity: (k+1+m+r) M_{r+1} + [(k+1)(1-x) - m x - r(2x-1)] M_r - r x(1-x) M_{r-1} = 0.
At x = y_j = (k+j)/(k+j+m), r = j-1:  (k+m+j) M_j = (j-1) y_j [M_{j-1} + (1-y_j) M_{j-2}].
Then the induction claim M_j(y_j) >= 0 for odd j, checked exactly."""
import random
from fractions import Fraction as Fr
from math import comb

def M(k, m, x, r):
    # exact: expand (x-y)^r and (1-y)^(m-1); int y^p dy = 1/(p+1)
    s = Fr(0)
    for i in range(r + 1):
        for l in range(m):
            p = k + (r - i) + l
            s += comb(r, i) * x ** i * (-1) ** (r - i) * comb(m - 1, l) * (-1) ** l * Fr(1, p + 1)
    return s

random.seed(1)
for _ in range(300):
    k, m, r = random.randint(0, 6), random.randint(1, 6), random.randint(1, 7)
    x = Fr(random.randint(0, 50), 50)
    lhs = (k + 1 + m + r) * M(k, m, x, r + 1) + ((k + 1) * (1 - x) - m * x - r * (2 * x - 1)) * M(k, m, x, r) \
        - r * x * (1 - x) * M(k, m, x, r - 1)
    assert lhs == 0, (k, m, r, x, lhs)
print('general recursion: 300 random exact checks OK')
for k in range(0, 8):
    for m in range(1, 8):
        for j in range(2, 12):
            y = Fr(k + j, k + j + m)
            assert (k + m + j) * M(k, m, y, j) == (j - 1) * y * (M(k, m, y, j - 1) + (1 - y) * M(k, m, y, j - 2))
        for j in range(1, 16, 2):
            assert M(k, m, Fr(k + j, k + j + m), j) >= 0
print('specialised recursion and M_j(y_j) >= 0 (odd j <= 15, k <= 7, m <= 7): OK')
