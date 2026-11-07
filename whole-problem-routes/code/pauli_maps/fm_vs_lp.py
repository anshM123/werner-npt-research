"""FM's GHZ necessary condition vs exact decomposability (LP over all Bell strings) for random non-trivial Pauli maps."""
import numpy as np, itertools, sys
from scipy.optimize import linprog
from pauli_general import bell, gam, fm_ok
n = int(sys.argv[1]); trials = int(sys.argv[2]); rng = np.random.default_rng(int(sys.argv[3]))
strings = list(itertools.product(range(4), repeat=n)); N = len(strings)
Gn = np.array([[np.prod([gam[J[s], K[s]] for s in range(n)]) for K in strings] for J in strings])
cnt = {}; mism = []
for _ in range(trials):
    l = rng.uniform(-1, 1, 3)
    f = bell(l)
    if f.min() >= 0 or (gam @ f).min() >= 0: continue
    F = np.array([np.prod(f[list(J)]) for J in strings])
    res = linprog(np.zeros(N), A_ub=Gn, b_ub=F, bounds=[(0, None)] * N, method='highs')
    lpok = res.status == 0; fm = fm_ok(l, n)
    cnt[(fm, lpok)] = cnt.get((fm, lpok), 0) + 1
    if fm != lpok: mism.append(np.round(l, 5))
print(f'n = {n}: (FM GHZ condition, decomposable) counts {cnt}; mismatches {mism[:5]}')
