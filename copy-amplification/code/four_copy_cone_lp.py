"""Four (n) copies from the three-copy theorem by convex decomposition.
Per site write a twirled Werner operator by its (P^s, P^a) eigenvalues normalised to (1, r):
   endpoint rho_{1/2}  ~ P^s + 3 P^a   (r = 3),
   separable Werner    ~ P^s + b P^a   with 0 <= b <= (1 + 1/d)/(1 - 1/d)   (d = 3: b <= 2),
   target rho_alpha    ~ P^s + rho P^a, rho = (1 + alpha)/(1 - alpha).
A product of per-site operators with at most m endpoint factors and separable others is undistillable whenever the
m-copy endpoint theorem holds (separable factors are absorbed into the local filters; d = 3: m <= 3 by the
deposited (unrefereed) three-copy certificates).  Mixtures of undistillable operators are undistillable.  Hence rho_alpha^{(x)n}
is n-copy undistillable if (1, rho, rho^2, ..., rho^n) lies in the cone of symmetrised generator moment vectors
   ( e_k(r_1..r_n) / C(n, k) )_{k = 0..n},   r_i in {3} u [0, bmax],   at most m of the r_i equal to 3.
(Necessity of this form: the products are diagonal in the parity basis Pi_A; symmetrising over site permutations
leaves the averages over |A| = k.)
We bisect on alpha for the largest feasible value.
usage: four_copy_cone_lp.py n m d
"""
import itertools
import sys
from math import comb

import numpy as np
from scipy.optimize import linprog

n, m, d = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
bmax = (1 + 1 / d) / (1 - 1 / d)
grid = list(np.linspace(0, bmax, 41)) + [3.0]
gens = []
for rs in itertools.combinations_with_replacement(grid, n):
    if sum(1 for r in rs if r == 3.0) > m:
        continue
    e = [sum(np.prod(c) for c in itertools.combinations(rs, k)) / comb(n, k) if k else 1.0 for k in range(n + 1)]
    gens.append(e)
G = np.array(gens).T          # (n+1) x #gens


def feasible(alpha):
    rho = (1 + alpha) / (1 - alpha)
    target = np.array([rho ** k for k in range(n + 1)])
    res = linprog(np.zeros(G.shape[1]), A_eq=G, b_eq=target, bounds=[(0, None)] * G.shape[1], method='highs')
    return res.status == 0


lo, hi = 1 / d, 0.5
assert feasible(lo)
for _ in range(40):
    mid = (lo + hi) / 2
    if feasible(mid):
        lo = mid
    else:
        hi = mid
print(f'n={n}, endpoint theorem up to m={m} copies, d={d}: rho_alpha^(x){n} undistillable for alpha <= {lo:.6f} '
      f'(NPT threshold 1/d = {1/d:.4f}; endpoint 1/2)')
