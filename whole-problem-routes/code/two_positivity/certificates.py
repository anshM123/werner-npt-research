"""Exact certificates (Proposition 7.11 of whole-problem-routes/README.md): for n copies and rational a,

    C_a^{(x)n} = sum_i lam_i  M_{p_i} (x) S_{i,1} (x) ... (x) S_{i,n-1}   +   sum_J r_J P_J ,     lam_i >= 0, r_J >= 0,

where M_p sits on copy p, every S is the O(3)-twirl of one of four product states (separable), and P_J runs over the
3^n products of P_Om, P_s0, P_an.  All operators involved are O(3)^n-invariant, so the identity is an identity between
coefficient vectors in the basis {P_J} (the P_J are orthogonal projectors, linearly independent); it is checked here
in exact rational arithmetic.  The weights lam_i are found by a floating-point LP (maximal slack) and then rounded to
rationals; nothing rests on the floating point: the rounded weights are re-checked exactly.

usage: python cert_ncopies.py n p q        (a = p/q)
"""
import sys, itertools
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog

n, p, q = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
a = Fr(p, q)
# coefficient vectors in the basis (P_Om, P_s0, P_an)
c = [a - 1, a - Fr(1, 4), a + Fr(1, 4)]                       # C_a
m = [Fr(-1, 2), Fr(1), Fr(1)]                                  # M = I - Omega Omega^*/2
# twirls of |e1 e1>, |w w> (w = (e1 + i e2)/sqrt2), |w conj(w)>, |e1 e2>: (t_Om / 1, t_s0 / 5, t_an / 3)
seps = {'e1e1': [Fr(1, 3), Fr(2, 15), Fr(0)], 'ww': [Fr(0), Fr(1, 5), Fr(0)],
        'wwbar': [Fr(1, 3), Fr(1, 30), Fr(1, 6)], 'e1e2': [Fr(0), Fr(1, 10), Fr(1, 6)]}


def kron(vs):
    out = [Fr(1)]
    for v in vs:
        out = [x * y for x in out for y in v]
    return out


target = kron([c] * n)
blocks, labels = [], []
for pos in range(n):
    for ss in itertools.product(sorted(seps), repeat=n - 1):
        vs = [seps[s] for s in ss]
        vs.insert(pos, m)
        blocks.append(kron(vs))
        labels.append((pos, ss))
B = np.array([[float(x) for x in b] for b in blocks]).T
t = np.array([float(x) for x in target])
# maximise eps subject to t - B lam >= eps, lam >= 0, sum lam <= 1e3
nb = B.shape[1]
cobj = np.zeros(nb + 1); cobj[-1] = -1
A_ub = np.hstack([B, np.ones((B.shape[0], 1))])
res = linprog(cobj, A_ub=A_ub, b_ub=t, bounds=[(0, None)] * nb + [(None, 1)], method='highs')
eps = -res.fun
print(f'n = {n}, a = {a} = {float(a):.6f}: floating-point max slack eps = {eps:.3e}')
if eps <= 0:
    print('no certificate at this a (FAIL)')
    sys.exit(1)
lam = [Fr(x).limit_denominator(10 ** 7) if x > 1e-12 else Fr(0) for x in res.x[:nb]]
# exact check
rem = list(target)
for l, b in zip(lam, blocks):
    if l:
        rem = [r - l * x for r, x in zip(rem, b)]
ok = all(l >= 0 for l in lam) and all(r >= 0 for r in rem)
used = [(labels[i], lam[i]) for i in range(nb) if lam[i]]
print(f'blocks used: {len(used)} of {nb}; min remainder coefficient = {float(min(rem)):.3e}')
for (pos, ss), l in used[:12]:
    print(f'   M on copy {pos + 1}, others {ss}: lam = {l}')
if len(used) > 12:
    print('   ...')
with open(f'certs/cert_n{n}_{p}_{q}.txt', 'w') as fh:
    fh.write(f'# n = {n}, a = {a}; lines: copy carrying M (1-based), separable factors on the other copies, lam\n')
    for (pos, ss), l in used:
        fh.write(f'{pos + 1} {" ".join(ss)} {l}\n')
print('exact check: all lam >= 0 and all remainder coefficients r_J >= 0:', ok)
print('PASS' if ok else 'FAIL')
