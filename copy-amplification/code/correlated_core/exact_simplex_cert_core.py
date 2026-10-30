"""[variant: env CORE_ALPHA = a0 sets the core rho_{a0}]
Exact (rational) Phase-I simplex for correlated-budget certificates at larger n.

Same generators and target as corr_budget_cert.py (factors E, O, X(b), Sig; at most m factors E; target
Binomial(n, p_alpha)).  Solves  G lam = target, sum lam = 1, lam >= 0  exactly with a dense rational tableau and
Bland's rule (no floating point at all), then writes corr_cert_d{d}_m{m}_a{A0.numerator}_{A0.denominator}_n{n}.json in the usual format.
usage: exact_simplex_cert.py d m n alpha_num alpha_den
"""
import json
import sys
import time
from fractions import Fraction as Fr
from math import comb

d, m, n = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
alpha = Fr(int(sys.argv[4]), int(sys.argv[5]))
import os
A0 = Fr(os.environ.get('CORE_ALPHA', '1/2'))   # core rho_{A0}, assumed m-copy undistillable
R0 = (1 + A0) / (1 - A0)
Ds, Da = Fr(d * (d + 1), 2), Fr(d * (d - 1), 2)
FAC = {'E': [Ds / (Ds + R0 * Da), R0 * Da / (Ds + R0 * Da)], 'O': [Fr(1), Fr(0)], 'X': [Fr(1, 2), Fr(1, 2)],
       'Sig': [Ds / d ** 2, Fr(0), Da / d ** 2]}


def conv(a, c):
    out = [Fr(0)] * (len(a) + len(c) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(c):
                out[i + j] += x * y
    return out


gens = []
for e in range(min(m, n) + 1):
    for s in range((n - e) // 2 + 1):
        rest = n - e - 2 * s
        for o in range(rest + 1):
            gens.append({'E': e, 'Sig': s, 'O': o, 'X': rest - o})
dists = []
for g in gens:
    dd = [Fr(1)]
    for name, c in g.items():
        for _ in range(c):
            dd = conv(dd, FAC[name])
    dists.append(dd)
p = Da * (1 + alpha) / (Ds * (1 - alpha) + Da * (1 + alpha))
target = [comb(n, k) * p ** k * (1 - p) ** (n - k) for k in range(n + 1)]
R = n + 2                                   # distribution rows + normalisation row
NG = len(gens)
# tableau rows: [A | I_R | b]; objective row: phase-I cost = sum of artificials
T = []
for r in range(R):
    row = [dists[j][r] if r <= n else Fr(1) for j in range(NG)]
    row += [Fr(1) if i == r else Fr(0) for i in range(R)]
    row.append(target[r] if r <= n else Fr(1))
    T.append(row)
basis = [NG + r for r in range(R)]
# reduced-cost row for minimising sum of artificials: z = -sum of constraint rows on the structural columns
obj = [-sum(T[r][j] for r in range(R)) for j in range(NG)] + [Fr(0)] * R + [-sum(T[r][-1] for r in range(R))]
t0 = time.time()
it = 0
while True:
    enter = next((j for j in range(NG + R) if obj[j] < 0), None)      # Bland: smallest index with negative cost
    if enter is None:
        break
    best, leave = None, None
    for r in range(R):
        a = T[r][enter]
        if a > 0:
            ratio = T[r][-1] / a
            if best is None or ratio < best or (ratio == best and basis[r] < basis[leave]):
                best, leave = ratio, r
    if leave is None:
        raise SystemExit('unbounded (cannot happen in phase I)')
    piv = T[leave][enter]
    T[leave] = [x / piv for x in T[leave]]
    for r in range(R):
        if r != leave and T[r][enter] != 0:
            f = T[r][enter]
            T[r] = [x - f * y for x, y in zip(T[r], T[leave])]
    f = obj[enter]
    obj = [x - f * y for x, y in zip(obj, T[leave])]
    basis[leave] = enter
    it += 1
    if it % 25 == 0:
        print(f'  pivot {it}: phase-I value {-obj[-1]} ~ {float(-obj[-1]):.3e} ({time.time() - t0:.0f}s)', flush=True)
phase1 = -obj[-1]
lam = [Fr(0)] * NG
for r, bv in enumerate(basis):
    if bv < NG:
        lam[bv] = T[r][-1]
feasible = phase1 == 0
ok = feasible and all(sum(dists[j][k] * lam[j] for j in range(NG)) == target[k] for k in range(n + 1)) \
    and sum(lam) == 1 and all(x >= 0 for x in lam)
used = [j for j in range(NG) if lam[j] != 0]
print(f'd={d} m={m} n={n} alpha={alpha} ({float(alpha):.6f}): phase-I optimum {float(phase1):.3e}; '
      f'{"CERTIFICATE VALID" if ok else "INFEASIBLE / FAILED"} with {len(used)} generators, {it} pivots, '
      f'{time.time() - t0:.0f}s', flush=True)
if ok:
    json.dump({'d': d, 'm': m, 'n': n, 'alpha': [alpha.numerator, alpha.denominator],
               'generators': [gens[j] for j in used], 'weights': [str(lam[j]) for j in used]},
              open(f'corr_cert_d{d}_m{m}_a{A0.numerator}_{A0.denominator}_n{n}.json', 'w'))
