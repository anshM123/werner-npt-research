"""n copies from <= m endpoint cores with named factors, including correlated separable blocks.

After the Werner twirl on every site and symmetrisation over sites, an operator on n sites is fixed by the
distribution of the number of antisymmetric sites; for a product of factors on disjoint sites that distribution is
the convolution of the factors' distributions.  Factors (per-factor distribution of #antisymmetric sites):
    E   (core, 1 site):       (D_s, 3 D_a)/(D_s + 3 D_a)
    O   (separable, 1 site):  (1, 0)
    X   (separable, 1 site):  (1/2, 1/2)                       X = X(b), the PPT-boundary Werner operator
    Sk  (separable, k sites): Alice holds GHZ_k on her k sites, Bob likewise (k = 2 is Sig); distribution computed
                              from Tr[rho_U^2] = 1/d for proper nonempty U (and 1 for U empty or everything).
Target: Binomial(n, p_alpha), p_alpha = D_a (1+alpha) / (D_s (1-alpha) + D_a (1+alpha)).
Feasibility LP over all multisets of factors (at most m factors E); bisection for the largest alpha.
usage: named_factor_ncopy.py d m kmax n1 n2 ...      (kmax = largest GHZ block used; 1 = no correlated blocks)
"""
import itertools
import sys
from fractions import Fraction as Fr
from math import comb

import numpy as np
from scipy.optimize import linprog

d, m, kmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
ns = [int(x) for x in sys.argv[4:]]
Ds, Da = d * (d + 1) / 2, d * (d - 1) / 2


def ghz_dist(k):
    """distribution of #antisymmetric sites for twirl(GHZ_k (x) GHZ_k) on k sites."""
    # p_T = 2^-k sum_U (-1)^{|U cap T|} tau_U,  tau_U = 1 (U empty or full), 1/d otherwise
    dist = np.zeros(k + 1)
    for t in range(k + 1):
        # weight of one fixed T with |T| = t, times C(k, t)
        T = set(range(t))
        tot = 0.0
        for r in range(k + 1):
            for U in itertools.combinations(range(k), r):
                tau = 1.0 if r in (0, k) else 1.0 / d
                tot += (-1) ** len(T & set(U)) * tau
        dist[t] = comb(k, t) * tot / 2 ** k
    return dist


F = {'E': np.array([Ds, 3 * Da]) / (Ds + 3 * Da), 'O': np.array([1.0, 0.0]), 'X': np.array([0.5, 0.5])}
for k in range(2, kmax + 1):
    F[f'S{k}'] = ghz_dist(k)
size = {'E': 1, 'O': 1, 'X': 1}
size.update({f'S{k}': k for k in range(2, kmax + 1)})


def generators(n):
    names = list(F)
    gens = []

    def rec(i, remaining, counts):
        if i == len(names):
            if remaining == 0 and counts.get('E', 0) <= m:
                dist = np.array([1.0])
                for nm, c in counts.items():
                    for _ in range(c):
                        dist = np.convolve(dist, F[nm])
                gens.append((dict(counts), dist))
            return
        nm = names[i]
        for c in range(remaining // size[nm] + 1):
            if nm == 'E' and c > m:
                break
            counts[nm] = c
            rec(i + 1, remaining - c * size[nm], counts)
        counts.pop(nm, None)
    rec(0, n, {})
    return gens


def target(n, alpha):
    p = Da * (1 + alpha) / (Ds * (1 - alpha) + Da * (1 + alpha))
    return np.array([comb(n, k) * p ** k * (1 - p) ** (n - k) for k in range(n + 1)])


def feasible(gens, n, alpha):
    G = np.array([g[1] for g in gens]).T
    res = linprog(np.zeros(G.shape[1]), A_eq=np.vstack([G, np.ones(G.shape[1])]),
                  b_eq=np.concatenate([target(n, alpha), [1.0]]), bounds=[(0, None)] * G.shape[1], method='highs')
    return res.status == 0


def budget_alpha(n):
    a0 = 0.5
    return (n * (1 + a0) + m * (d * a0 - 1)) / (n * d * (1 + a0) - m * (d * a0 - 1))


for n in ns:
    gens = generators(n)
    lo, hi = 1.0 / d, 0.5
    for _ in range(45):
        mid = (lo + hi) / 2
        if feasible(gens, n, mid):
            lo = mid
        else:
            hi = mid
    b = budget_alpha(n)
    print(f'd={d} m={m} kmax={kmax} n={n}: threshold {lo:.6f}  (budget {b:.6f}; n*(alpha-1/d): '
          f'{n * (lo - 1 / d):.4f} vs {n * (b - 1 / d):.4f}); {len(gens)} generator types', flush=True)
