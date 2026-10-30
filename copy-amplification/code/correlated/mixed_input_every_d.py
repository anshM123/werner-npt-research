"""Four copies in every dimension from MIXED inputs: two-copy endpoint theorem and MacCormack-Bao three-copy bound.

Cores: E = X(3) (endpoint) and G = X(117/43) (alpha_0 = 37/80).  Allowed core patterns (each one-copy undistillable):
    E^e with e <= 2                       (two-copy theorem),
    E (x) G                               (G is a mixture of E and the separable X(b), so two-copy theorem + absorption),
    G^g with g <= 3                       (MacCormack-Bao),
times separable factors O, X(b) and the correlated pair Sig.  Threshold by LP bisection (n = 4), compared with the
single-input thresholds.
usage: mixed_input_every_d.py d1 d2 ...
"""
import sys
from fractions import Fraction as Fr
from math import comb

import numpy as np
from scipy.optimize import linprog

n = 4
r_E, r_G = 3.0, 117 / 43


def dists(d, patterns):
    Ds, Da = d * (d + 1) / 2, d * (d - 1) / 2
    F = {'E': np.array([Ds, r_E * Da]) / (Ds + r_E * Da), 'G': np.array([Ds, r_G * Da]) / (Ds + r_G * Da),
         'O': np.array([1.0, 0.0]), 'X': np.array([0.5, 0.5]), 'Sig': np.array([Ds, 0.0, Da]) / d ** 2}
    out = []
    for (e, g) in patterns:
        rest_sites = n - e - g
        for s in range(rest_sites // 2 + 1):
            rest = rest_sites - 2 * s
            for o in range(rest + 1):
                dist = np.array([1.0])
                for name, cnt in (('E', e), ('G', g), ('Sig', s), ('O', o), ('X', rest - o)):
                    for _ in range(cnt):
                        dist = np.convolve(dist, F[name])
                out.append(dist)
    return np.array(out).T


def threshold(d, patterns):
    Ds, Da = d * (d + 1) / 2, d * (d - 1) / 2
    G = dists(d, patterns)
    lo, hi = 1.0 / d, 0.5
    for _ in range(45):
        a = (lo + hi) / 2
        p = Da * (1 + a) / (Ds * (1 - a) + Da * (1 + a))
        t = np.array([comb(n, k) * p ** k * (1 - p) ** (n - k) for k in range(n + 1)])
        res = linprog(np.zeros(G.shape[1]), A_eq=np.vstack([G, np.ones(G.shape[1])]), b_eq=np.concatenate([t, [1.0]]),
                      bounds=[(0, None)] * G.shape[1], method='highs')
        if res.status == 0:
            lo = a
        else:
            hi = a
    return lo


two = [(e, 0) for e in range(3)]
three_G = [(0, g) for g in range(4)]
mixed = two + three_G + [(1, 1)]
for d in [int(x) for x in sys.argv[1:]]:
    t2, t3, tm = threshold(d, two), threshold(d, three_G), threshold(d, mixed)
    print(f'd={d}: two-copy only {t2:.5f} | MB only {t3:.5f} | mixed {tm:.5f} (gain over best {tm - max(t2, t3):+.5f})',
          flush=True)
