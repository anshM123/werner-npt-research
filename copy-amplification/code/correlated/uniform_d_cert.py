"""Dimension-uniform four-copy certificate from the two-copy theorem: symbolic in d.

Fix alpha (rational).  For d large, the four-copy target and the factor distributions (G = endpoint core,
O, X(b), Sig_d) are rational functions of d.  Take a support set of generators (found numerically at large d), solve
the linear system for the weights SYMBOLICALLY in d, and prove that every weight is >= 0 for all d >= D0 (sign
analysis of numerator and denominator polynomials via exact real-root isolation).  Together with an exact certificate
at d = D0 (which covers every d <= D0, because compressing a Werner state to a subspace gives the Werner state with
the same alpha), this proves four-copy undistillability for alpha <= alpha in EVERY dimension.
usage: uniform_d_cert.py alpha_num alpha_den D0 [m r0_num r0_den [n]]
"""
import sys
from fractions import Fraction as Fr
from math import comb

import numpy as np
import sympy as sp
from scipy.optimize import linprog

an, ad, D0 = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
M_ARG = int(sys.argv[4]) if len(sys.argv) > 4 else 2
R0 = sp.Rational(int(sys.argv[5]), int(sys.argv[6])) if len(sys.argv) > 6 else sp.Integer(3)
alpha = sp.Rational(an, ad)
m, n = M_ARG, (int(sys.argv[7]) if len(sys.argv) > 7 else 4)
d = sp.symbols('d', positive=True)


def factor_dists(dv):
    Ds, Da = dv * (dv + 1) / 2, dv * (dv - 1) / 2
    return {'G': [Ds / (Ds + R0 * Da), R0 * Da / (Ds + R0 * Da)], 'O': [1, 0], 'X': [sp.Rational(1, 2), sp.Rational(1, 2)],
            'Sig': [Ds / dv ** 2, 0, Da / dv ** 2]}


gens = []
for e in range(m + 1):
    for s in range((n - e) // 2 + 1):
        rest = n - e - 2 * s
        for o in range(rest + 1):
            gens.append({'G': e, 'Sig': s, 'O': o, 'X': rest - o})


def dist(g, F):
    out = [1]
    for name, c in g.items():
        for _ in range(c):
            new = [0] * (len(out) + len(F[name]) - 1)
            for i, x in enumerate(out):
                for j, y in enumerate(F[name]):
                    new[i + j] += x * y
            out = new
    return out


def target(dv):
    Ds, Da = dv * (dv + 1) / 2, dv * (dv - 1) / 2
    p = Da * (1 + alpha) / (Ds * (1 - alpha) + Da * (1 + alpha))
    return [comb(n, k) * p ** k * (1 - p) ** (n - k) for k in range(n + 1)]


# numeric support at a large d
dnum = 10 ** 6
Fn = factor_dists(sp.Integer(dnum))
G = np.array([[float(x) for x in dist(g, Fn)] for g in gens]).T
t = np.array([float(x) for x in target(sp.Integer(dnum))])
res = linprog(np.zeros(G.shape[1]), A_eq=np.vstack([G, np.ones(G.shape[1])]), b_eq=np.concatenate([t, [1.0]]),
              bounds=[(0, None)] * G.shape[1], method='highs-ds')
assert res.status == 0, 'infeasible at large d'
support = [i for i in range(len(gens)) if res.x[i] > 1e-12]
print('support at d = 1e6:', [gens[i] for i in support])
# symbolic solve on the support (n+1 distribution equations + normalisation; one is redundant)
Fs = factor_dists(d)
lam = sp.symbols(f'l0:{len(support)}')
eqs = []
tg = target(d)
for k in range(n + 1):
    eqs.append(sp.Eq(sum(lam[j] * dist(gens[support[j]], Fs)[k] for j in range(len(support))), tg[k]))
eqs.append(sp.Eq(sum(lam), 1))
sol = sp.solve(eqs, lam, dict=True)
assert sol, 'no symbolic solution on this support'
sol = sol[0]
ok_all = True
for j in range(len(support)):
    w = sp.factor(sp.together(sol[lam[j]]))
    num, den = sp.fraction(w)
    # sign on [D0, oo): check no real roots of num, den in [D0, oo) and sign at a sample point
    bad = False
    for poly in (num, den):
        P = sp.Poly(sp.expand(poly), d)
        roots = [r for r in P.real_roots() if r >= D0]
        if roots:
            bad = True
    val = w.subs(d, D0 + 1)
    nonneg = (not bad) and val >= 0
    ok_all &= nonneg
    print(f'  weight {gens[support[j]]}: sign-constant on [{D0}, oo): {not bad}; value at d={D0 + 1}: {float(val):.4e}')
print(f'ALL WEIGHTS >= 0 FOR EVERY d >= {D0}: {ok_all}  (alpha = {alpha} = {float(alpha)})')
with open(f'uniform_d_symbolic_weights_{an}_{ad}_D{D0}' + ('' if R0 == 3 else f'_m{m}_core{R0.p}_{R0.q}') + ('' if n == 4 else f'_n{n}') + '.txt', 'w') as fh:
    fh.write(f'alpha = {alpha}; every weight >= 0 for all d >= {D0}: {ok_all}\n')
    fh.write('generator (factor counts): weight as a rational function of d\n')
    for j in range(len(support)):
        fh.write(f'{gens[support[j]]}: {sp.factor(sp.together(sol[lam[j]]))}\n')
