"""Correlated-factor windows in every dimension (four copies by default), with exact rational certificates.

Core G = X(r0) = P^s + r0 P^a, assumed m-copy undistillable (r0 = 3 for the endpoint; r0 = 117/43 for alpha_0 = 37/80,
MacCormack-Bao's three-copy bound in every d).  Separable factors O = P^s, X(b) (b = (d+1)/(d-1)) and the correlated
pair Sig_d (distribution (D_s, 0, D_a)/d^2).  Threshold by float-LP bisection, then an exact rational Phase-I simplex
certifies a rational alpha just below it.  Also prints the product budget and the mean-bound ceiling.
usage: corr_every_d.py m r0_num r0_den n d1 d2 ...
"""
import sys
from fractions import Fraction as Fr
from math import comb

import numpy as np
from scipy.optimize import linprog

m = int(sys.argv[1]); r0 = Fr(int(sys.argv[2]), int(sys.argv[3])); n = int(sys.argv[4])
ds = [int(x) for x in sys.argv[5:]]


def conv(a, c):
    out = [Fr(0)] * (len(a) + len(c) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(c):
                out[i + j] += x * y
    return out


def setup(d):
    Ds, Da = Fr(d * (d + 1), 2), Fr(d * (d - 1), 2)
    FAC = {'G': [Ds / (Ds + r0 * Da), r0 * Da / (Ds + r0 * Da)], 'O': [Fr(1), Fr(0)], 'X': [Fr(1, 2), Fr(1, 2)],
           'Sig': [Ds / d ** 2, Fr(0), Da / d ** 2]}
    gens = []
    for e in range(min(m, n) + 1):
        for s in range((n - e) // 2 + 1):
            rest = n - e - 2 * s
            for o in range(rest + 1):
                gens.append({'G': e, 'Sig': s, 'O': o, 'X': rest - o})
    dists = []
    for g in gens:
        dd = [Fr(1)]
        for name, c in g.items():
            for _ in range(c):
                dd = conv(dd, FAC[name])
        dists.append(dd)
    return Ds, Da, gens, dists


def target(Ds, Da, alpha):
    p = Da * (1 + alpha) / (Ds * (1 - alpha) + Da * (1 + alpha))
    return [comb(n, k) * p ** k * (1 - p) ** (n - k) for k in range(n + 1)]


def exact_feasible(dists, tgt):
    R = n + 2; NG = len(dists)
    T = []
    for r in range(R):
        row = [dists[j][r] if r <= n else Fr(1) for j in range(NG)] + [Fr(1) if i == r else Fr(0) for i in range(R)]
        row.append(tgt[r] if r <= n else Fr(1))
        T.append(row)
    basis = [NG + r for r in range(R)]
    obj = [-sum(T[r][j] for r in range(R)) for j in range(NG)] + [Fr(0)] * R + [-sum(T[r][-1] for r in range(R))]
    while True:
        enter = next((j for j in range(NG + R) if obj[j] < 0), None)
        if enter is None:
            break
        best, leave = None, None
        for r in range(R):
            a = T[r][enter]
            if a > 0:
                ratio = T[r][-1] / a
                if best is None or ratio < best or (ratio == best and basis[r] < basis[leave]):
                    best, leave = ratio, r
        piv = T[leave][enter]
        T[leave] = [x / piv for x in T[leave]]
        for r in range(R):
            if r != leave and T[r][enter] != 0:
                f = T[r][enter]
                T[r] = [x - f * y for x, y in zip(T[r], T[leave])]
        f = obj[enter]
        obj = [x - f * y for x, y in zip(obj, T[leave])]
        basis[leave] = enter
    lam = [Fr(0)] * NG
    for r, bv in enumerate(basis):
        if bv < NG:
            lam[bv] = T[r][-1]
    ok = (-obj[-1] == 0) and all(sum(dists[j][k] * lam[j] for j in range(NG)) == tgt[k] for k in range(n + 1)) \
        and sum(lam) == 1 and all(x >= 0 for x in lam)
    return ok, lam


for d in ds:
    Ds, Da, gens, dists = setup(d)
    G = np.array([[float(x) for x in dd] for dd in dists]).T
    lo, hi = 1.0 / d, 0.5
    for _ in range(45):
        mid = (lo + hi) / 2
        t = np.array([float(x) for x in target(Ds, Da, Fr(mid))])
        res = linprog(np.zeros(G.shape[1]), A_eq=np.vstack([G, np.ones(G.shape[1])]), b_eq=np.concatenate([t, [1.0]]),
                      bounds=[(0, None)] * G.shape[1], method='highs')
        if res.status == 0:
            lo = mid
        else:
            hi = mid
    a_cert = Fr(int(lo * 10 ** 5) - 1, 10 ** 5)
    ok, lam = exact_feasible(dists, target(Ds, Da, a_cert))
    if ok:
        import json
        used = [j for j in range(len(gens)) if lam[j] != 0]
        json.dump({'d': d, 'm': m, 'n': n, 'core_r0': [r0.numerator, r0.denominator],
                   'alpha': [a_cert.numerator, a_cert.denominator],
                   'generators': [{('E' if k == 'G' else k): v for k, v in gens[j].items()} for j in used],
                   'weights': [str(lam[j]) for j in used]}, open(f'corr_cert_d{d}_m{m}_n{n}' + ('' if r0 == 3 else f'_core{r0.numerator}_{r0.denominator}') + '.json', 'w'))
    b = Fr(d + 1, d - 1)
    a0 = (r0 - 1) / (r0 + 1)
    budget = (n * (1 + a0) + m * (d * a0 - 1)) / (n * d * (1 + a0) - m * (d * a0 - 1))
    R_ = Fr(m, n) * (r0 - b) / (r0 + b)
    rho = b * (1 + R_) / (1 - R_)
    ceiling = (rho - 1) / (rho + 1)
    print(f'd={d}: threshold {lo:.6f}; exact certificate at {a_cert} ({float(a_cert):.5f}): {"VALID" if ok else "FAILED"}; '
          f'budget {float(budget):.6f}; ceiling {float(ceiling):.6f}; gain/(ceiling-budget) = '
          f'{(lo - float(budget)) / float(ceiling - budget):.2f}', flush=True)
