"""Independent exact check of the correlated-budget certificates with a NON-endpoint core (files
corr_cert_d{d}_m{m}_a{p}_{q}_n{n}.json).  A certificate is valid iff: every generator uses at most m cores E and fills
n sites, all weights are >= 0, the weights sum to 1, and the mixture of generator distributions equals the
Binomial(n, p_alpha) target exactly.  Core E = rho_{a0} with ratio (1+a0)/(1-a0).  Also checks the closed-form
all-n windows obtained from the n = 24 certificates through the budget theorem."""
import glob
import json
import re
from fractions import Fraction as Fr
from math import comb


def check(path):
    c = json.load(open(path))
    d, m, n = c['d'], c['m'], c['n']
    p, q = map(int, re.search(r'_a(\d+)_(\d+)_n', path).groups()); a0 = Fr(p, q)
    alpha = Fr(*c['alpha']); w = [Fr(x) for x in c['weights']]
    Ds, Da = Fr(d * (d + 1), 2), Fr(d * (d - 1), 2); R0 = (1 + a0) / (1 - a0)
    FAC = {'E': [Ds / (Ds + R0 * Da), R0 * Da / (Ds + R0 * Da)], 'O': [Fr(1), Fr(0)], 'X': [Fr(1, 2), Fr(1, 2)],
           'Sig': [Ds / d ** 2, Fr(0), Da / d ** 2]}
    SIZE = {'E': 1, 'O': 1, 'X': 1, 'Sig': 2}
    tot = [Fr(0)] * (n + 1)
    for g, wt in zip(c['generators'], w):
        if g.get('E', 0) > m or sum(SIZE[k] * v for k, v in g.items()) != n or wt < 0:
            return d, m, n, alpha, a0, False
        dist = [Fr(1)]
        for k, v in g.items():
            for _ in range(v):
                dist = [sum(dist[i] * FAC[k][j - i] for i in range(len(dist)) if 0 <= j - i < len(FAC[k]))
                        for j in range(len(dist) + len(FAC[k]) - 1)]
        for k_ in range(n + 1):
            tot[k_] += wt * dist[k_]
    pa = Da * (1 + alpha) / (Ds * (1 - alpha) + Da * (1 + alpha))
    ok = sum(w) == 1 and all(tot[k] == comb(n, k) * pa ** k * (1 - pa) ** (n - k) for k in range(n + 1))
    return d, m, n, alpha, a0, ok


rows = sorted(check(f) for f in glob.glob('corr_cert_d*_m*_a*_*_n*.json'))
for d, m, n, alpha, a0, ok in rows:
    print(f'd={d} core a0={a0} m={m} n={n:2d}: alpha={alpha} ({float(alpha):.5f})  {"VALID" if ok else "INVALID"}')
# closed forms for n >= 24 via the budget theorem with the n = 24 certificate as input
for d, a24, num, den in ((3, Fr(226, 625), (851, 1272), (2553, -1272)), (4, Fr(2913, 10000), (12913, 39648), (51652, -39648))):
    kp, k24 = Fr(d - 1, d + 1), (1 - a24) / (1 + a24)
    good = True
    for n in range(24, 400):
        kap = Fr(24, n) * k24 + (1 - Fr(24, n)) * kp
        good &= (1 - kap) / (1 + kap) == Fr(num[0] * n + num[1], den[0] * n + den[1])
    print(f'd={d}: all n >= 24 window ({num[0]} n + {num[1]})/({den[0]} n {den[1]:+d}) matches budget formula: {good};'
          f' n*(alpha_n - 1/d) -> {float(Fr(d * num[1] - den[1], d * den[0]) * d) / d:.4f}')
