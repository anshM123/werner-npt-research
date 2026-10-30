"""Standalone verifier for correlated-budget certificates corr_cert_d{d}_m{m}_n{n}.json (exact, pure Python Fractions).

Claim: rho_alpha^{(x)n} (Werner, rho ~ I - alpha F on each of n sites, local dimension d) is a convex combination of
Werner-twirled, site-symmetrised products of the named factors
    E = P^s + 3 P^a (endpoint, at most m of them), O = P^s, X = P^s + b P^a (b = (d+1)/(d-1)),
    Sig = |Phi><Phi|_{a_i a_j} (x) |Phi><Phi|_{b_i b_j}  (Alice and Bob each hold a maximally entangled pair; separable).
Such operators are fixed by the distribution of the number of antisymmetric sites; for products it is the convolution
of the factors' distributions, which this script recomputes from scratch:
    E: (D_s, 3 D_a)/(D_s + 3 D_a),  O: (1, 0),  X: (D_s, b D_a)/(D_s + b D_a) = (1/2, 1/2),
    Sig: <Phi_A Phi_B| P^x_i (x) P^y_j |Phi_A Phi_B> = Tr(P^x P^y)/d^2  ->  (D_s, 0, D_a)/d^2.
With --dense it also checks these per-factor weights against explicit matrices (small d).
Every generator is one-copy undistillable given the m-copy endpoint theorem (absorb separable factors), so a valid
certificate proves that rho_alpha is n-copy undistillable.
usage: verify_corr_cert.py cert.json [--dense]
"""
import itertools
import json
import sys
from fractions import Fraction as Fr
from math import comb

import numpy as np

cert = json.load(open(sys.argv[1]))
d, m, n = cert['d'], cert['m'], cert['n']
alpha = Fr(*cert['alpha'])
Ds, Da = Fr(d * (d + 1), 2), Fr(d * (d - 1), 2)
b = Fr(d + 1, d - 1)
r0 = Fr(*cert['core_r0']) if 'core_r0' in cert else Fr(3)          # core E = P^s + r0 P^a (r0 = 3: endpoint)
FAC = {'E': [Ds / (Ds + r0 * Da), r0 * Da / (Ds + r0 * Da)], 'O': [Fr(1), Fr(0)],
       'X': [Ds / (Ds + b * Da), b * Da / (Ds + b * Da)], 'Sig': [Ds / d ** 2, Fr(0), Da / d ** 2]}
SIZE = {'E': 1, 'O': 1, 'X': 1, 'Sig': 2}

if '--dense' in sys.argv:
    F = np.zeros((d * d, d * d))
    for i in range(d):
        for j in range(d):
            F[i * d + j, j * d + i] = 1
    Ps, Pa = (np.eye(d * d) + F) / 2, (np.eye(d * d) - F) / 2
    for name, op in (('E', Ps + float(r0) * Pa), ('O', Ps), ('X', Ps + float(b) * Pa)):
        w = np.array([np.trace(Ps @ op), np.trace(Pa @ op)]); w /= w.sum()
        assert np.allclose(w, [float(x) for x in FAC[name]]), name
    phi = np.eye(d).reshape(-1) / np.sqrt(d)                     # Phi on (a_i a_j) and on (b_i b_j)
    psi = np.kron(phi, phi).reshape(d, d, d, d)                  # indices a_i a_j b_i b_j
    psi = psi.transpose(0, 2, 1, 3).reshape(-1)                  # -> (a_i b_i)(a_j b_j)
    got = []
    for x, y in ((Ps, Ps), (Ps, Pa), (Pa, Ps), (Pa, Pa)):
        got.append(psi @ np.kron(x, y) @ psi)
    want = [float(Ds / d ** 2), 0.0, 0.0, float(Da / d ** 2)]
    assert np.allclose(got, want), (got, want)
    print('dense check of factor weights: OK (Sig weights ss, sa, as, aa =', np.round(got, 6), ')')


def conv(a, c):
    out = [Fr(0)] * (len(a) + len(c) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(c):
            out[i + j] += x * y
    return out


p = Da * (1 + alpha) / (Ds * (1 - alpha) + Da * (1 + alpha))
target = [comb(n, k) * p ** k * (1 - p) ** (n - k) for k in range(n + 1)]
total = [Fr(0)] * (n + 1)
lam = [Fr(x) for x in cert['weights']]
for g, w in zip(cert['generators'], lam):
    assert sum(SIZE[k] * v for k, v in g.items()) == n and g.get('E', 0) <= m, g
    dist = [Fr(1)]
    for name, c in g.items():
        for _ in range(c):
            dist = conv(dist, FAC[name])
    for k in range(n + 1):
        total[k] += w * dist[k]
ok = total == target and sum(lam) == 1 and all(x >= 0 for x in lam)
print(f'd={d} m={m} n={n} core r0={r0} alpha={alpha} ({float(alpha):.6f}): {len(lam)} generators, weights >= 0: '
      f'{all(x >= 0 for x in lam)}, sum = 1: {sum(lam) == 1}, identity exact: {total == target}')
print('CERTIFICATE VALID:', ok)
