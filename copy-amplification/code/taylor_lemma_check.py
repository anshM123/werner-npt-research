"""Lemma check: for t = n/(n-m), (1+tw)^n = P_{m-1}(w) + w^m R(w) with R having nonnegative coefficients in the
basis (1+w)^k.  Exact rational arithmetic; also scan t slightly below the bound."""
from fractions import Fraction as Fr
from math import comb
import sys

def coeffs_R_in_r(n, m, t):
    # R(w) = sum_{j>=m} C(n,j) t^j w^{j-m};  substitute w = r - 1 and collect powers of r
    out = [Fr(0)] * (n - m + 1)
    for j in range(m, n + 1):
        a = comb(n, j) * t ** j
        p = j - m                       # w^p = (r-1)^p
        for k in range(p + 1):
            out[k] += a * comb(p, k) * (-1) ** (p - k)
    return out

bad = []
for n in range(2, 31):
    for m in range(1, n):
        t = Fr(n, n - m)
        c = coeffs_R_in_r(n, m, t)
        if min(c) < 0:
            bad.append((n, m, min(c)))
print('cases with a negative coefficient at t = n/(n-m):', bad[:10], f'(total {len(bad)})')
