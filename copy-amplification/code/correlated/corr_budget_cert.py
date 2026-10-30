"""Exact certificates for n-copy windows with correlated separable pairs (named factors), any d, m, n.

Factors and their distributions of the number of antisymmetric sites (after the Werner twirl on every site):
    E   = X(3)                     (core)      : (D_s, 3 D_a) / (D_s + 3 D_a)
    O   = P^s                      (separable) : (1, 0)
    X   = X(b), b = (d+1)/(d-1)    (separable) : (1/2, 1/2)
    Sig = twirl(Phi_A (x) Phi_B)   (separable, two sites; Alice and Bob each hold a maximally entangled pair)
                                                 : (D_s, 0, D_a) / d^2
The twirled, site-symmetrised operator of a product of factors on disjoint sites is fixed by the convolution of these
distributions.  A generator is a multiset of factors filling n sites with at most m factors E.  If the target
Binomial(n, p_alpha), p_alpha = D_a(1+alpha)/(D_s(1-alpha) + D_a(1+alpha)), is a convex combination of generator
distributions, then rho_alpha^{(x)n} is a convex combination of local-unitary images of generators, each of which is
one-copy undistillable given the m-copy endpoint theorem (separable factors are absorbed).  Hence rho_alpha is
n-copy undistillable.
Finds the largest alpha by LP bisection, then certifies a rational alpha just below it EXACTLY (Fractions) and writes
the certificate to corr_cert_d{d}_m{m}_n{n}.json.
usage: corr_budget_cert.py d m n1 n2 ...
"""
import json
import sys
from fractions import Fraction as Fr
from math import comb

import numpy as np
from scipy.optimize import linprog

d, m = int(sys.argv[1]), int(sys.argv[2])
import os
MARGIN = float(os.environ.get("MARGIN", "0"))
ns = [int(x) for x in sys.argv[3:]]
Ds, Da = Fr(d * (d + 1), 2), Fr(d * (d - 1), 2)
FAC = {'E': [Ds / (Ds + 3 * Da), 3 * Da / (Ds + 3 * Da)], 'O': [Fr(1), Fr(0)], 'X': [Fr(1, 2), Fr(1, 2)],
       'Sig': [Ds / d ** 2, Fr(0), Da / d ** 2]}
SIZE = {'E': 1, 'O': 1, 'X': 1, 'Sig': 2}


def conv(a, b):
    out = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def generators(n):
    gens = []
    for e in range(min(m, n) + 1):
        for s in range((n - e) // 2 + 1):
            rest = n - e - 2 * s
            for o in range(rest + 1):
                x = rest - o
                gens.append({'E': e, 'Sig': s, 'O': o, 'X': x})
    return gens


def dist(g):
    out = [Fr(1)]
    for name, c in g.items():
        for _ in range(c):
            out = conv(out, FAC[name])
    return out


def target(n, alpha):
    p = Da * (1 + alpha) / (Ds * (1 - alpha) + Da * (1 + alpha))
    return [comb(n, k) * p ** k * (1 - p) ** (n - k) for k in range(n + 1)]


def lp(G, t, method='highs'):
    return linprog(np.zeros(G.shape[1]), A_eq=np.vstack([G, np.ones(G.shape[1])]), b_eq=np.concatenate([t, [1.0]]),
                   bounds=[(0, None)] * G.shape[1], method=method)


def budget(n):
    a0 = Fr(1, 2)
    return (n * (1 + a0) + m * (d * a0 - 1)) / (n * d * (1 + a0) - m * (d * a0 - 1))


for n in ns:
    gens = generators(n)
    D = [dist(g) for g in gens]
    G = np.array([[float(x) for x in Dk] for Dk in D]).T
    lo, hi = 1.0 / d, 0.5
    for _ in range(50):
        mid = (lo + hi) / 2
        if lp(G, np.array([float(x) for x in target(n, Fr(mid))])).status == 0:
            lo = mid
        else:
            hi = mid
    # certify at a rational alpha slightly below the threshold
    alpha = Fr(int((lo - MARGIN) * 10 ** 5) - 1, 10 ** 5)
    tF = target(n, alpha)
    res = lp(G, np.array([float(x) for x in tF]), method='highs-ipm')
    xf = res.x
    used = [i for i in range(len(gens)) if xf[i] > 1e-11]
    k = len(used)
    rows = [[D[j][r] for j in used] + [tF[r]] for r in range(n + 1)] + [[Fr(1)] * k + [Fr(1)]]
    # exact solve: Gaussian elimination; free (non-pivot) variables are fixed at their rationalised float values
    M = [r[:] for r in rows]; piv = []; r0 = 0
    for c in range(k):
        p_ = next((r for r in range(r0, len(M)) if M[r][c] != 0), None)
        if p_ is None:
            continue
        M[r0], M[p_] = M[p_], M[r0]
        for r in range(len(M)):
            if r != r0 and M[r][c] != 0:
                f = M[r][c] / M[r0][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[r0])]
        piv.append(c); r0 += 1
    free = [c for c in range(k) if c not in piv]
    lam = [Fr(0)] * k
    for c in free:
        lam[c] = Fr(float(xf[used[c]])).limit_denominator(10 ** 12)
    for i, c in enumerate(piv):
        lam[c] = (M[i][k] - sum(M[i][f] * lam[f] for f in free)) / M[i][c]
    ok = all(sum(rows[r][j] * lam[j] for j in range(k)) == rows[r][k] for r in range(len(rows)))
    valid = ok and all(x >= 0 for x in lam)
    b = budget(n)
    print(f'd={d} m={m} n={n}: threshold {lo:.6f} (budget {float(b):.6f}); exact certificate at alpha={alpha} '
          f'({float(alpha):.5f}): {"VALID" if valid else "FAILED"} with {k} generators', flush=True)
    json.dump({'d': d, 'm': m, 'n': n, 'alpha': [alpha.numerator, alpha.denominator],
               'generators': [gens[used[j]] for j in range(k)], 'weights': [str(x) for x in lam]},
              open(f'corr_cert_d{d}_m{m}_n{n}.json', 'w'))
