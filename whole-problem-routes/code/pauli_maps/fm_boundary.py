"""Near-boundary test: for random directions, bisect the FM GHZ condition along the ray from 0 (where Y is CP) to the
cube boundary, then test decomposability just inside (0.999 x) and just outside (1.001 x) the FM boundary."""
import numpy as np, itertools, sys
from scipy.optimize import linprog
from pauli_general import bell, gam, fm_ok
n = int(sys.argv[1]); trials = int(sys.argv[2]); rng = np.random.default_rng(int(sys.argv[3]))
strings = list(itertools.product(range(4), repeat=n)); N = len(strings)
Gn = np.array([[np.prod([gam[J[s], K[s]] for s in range(n)]) for K in strings] for J in strings])
def dec(l):
    f = bell(l); F = np.array([np.prod(f[list(J)]) for J in strings])
    return linprog(np.zeros(N), A_ub=Gn, b_ub=F, bounds=[(0, None)] * N, method='highs').status == 0
res = []
for _ in range(trials):
    d = rng.standard_normal(3); d /= np.abs(d).max()          # ray to the cube boundary
    lo, hi = 0.0, 1.0
    if fm_ok(d, n): continue                                    # whole ray fine: skip
    for _ in range(40):
        mid = (lo + hi) / 2
        if fm_ok(mid * d, n): lo = mid
        else: hi = mid
    inside, outside = dec(0.999 * lo * d), dec(1.001 * lo * d)
    trivial = (bell(0.999 * lo * d).min() >= 0) or ((gam @ bell(0.999 * lo * d)).min() >= 0)
    res.append((inside, outside, trivial))
ok = sum(1 for r in res if r[0] and not r[1])
print(f'n = {n}: rays {len(res)}; decomposable just inside and not just outside the FM boundary: {ok}; '
      f'others: {[r for r in res if not (r[0] and not r[1])][:5]}')
