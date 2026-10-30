"""d -> infinity limit of the correlated-factor LP (two-copy cores), to locate the limit of the method.

At d = infinity the factor distributions of the number of antisymmetric sites are
    G = endpoint core E = P^s + 3 P^a: Bern(3/4)  (at most m per term),  O = P^s: delta_0,  X(b) -> Bern(1/2),
    Sig (two sites): (1/2) delta_0 + (1/2) delta_2,
and rho_alpha^{(x)n} has Bin(n, (1+alpha)/2).  The ceiling (Prop. ceiling) reads n alpha <= m/2 here.
Bisection on alpha of the LP feasibility "target in the convex hull of the generator distributions".
usage: dinf_limit_lp.py m n1 n2 ...
"""
import sys
from math import comb

import numpy as np
from scipy.optimize import linprog

m = int(sys.argv[1])
FAC = {'G': [0.25, 0.75], 'O': [1.0, 0.0], 'X': [0.5, 0.5], 'Sig': [0.5, 0.0, 0.5]}


def gens(n):
    out = []
    for e in range(m + 1):
        for s in range((n - e) // 2 + 1):
            rest = n - e - 2 * s
            for o in range(rest + 1):
                dist = np.array([1.0])
                for name, c in (('G', e), ('Sig', s), ('O', o), ('X', rest - o)):
                    for _ in range(c):
                        dist = np.convolve(dist, FAC[name])
                out.append(dist)
    return np.array(out).T


def feasible(n, a, M):
    p = (1 + a) / 2
    t = np.array([comb(n, k) * p ** k * (1 - p) ** (n - k) for k in range(n + 1)])
    # scale rows to make tiny tail probabilities count
    w = 1 / np.maximum(t, 1e-300)
    res = linprog(np.zeros(M.shape[1]), A_eq=M * w[:, None], b_eq=t * w, bounds=(0, None), method='highs')
    return res.status == 0


for n in map(int, sys.argv[2:]):
    M = gens(n)
    lo, hi = 0.0, min(1.0, m / (2 * n) + 1e-9)
    for _ in range(40):
        mid = (lo + hi) / 2
        if feasible(n, mid, M):
            lo = mid
        else:
            hi = mid
    print(f'n={n:3d}: alpha_n(d=inf) ~ {lo:.6f}   n*alpha = {n * lo:.5f}   ceiling n*alpha <= {m / 2}   '
          f'dimension-free constant n*alpha/(1+alpha) = {n * lo / (1 + lo):.5f}', flush=True)
