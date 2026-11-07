"""Rational certificate from a sym_sdp.py solution (diagram coordinates xA, xB).

Rounds to integers over 2^log2L, then symmetrises EXACTLY over permutations of the copies and over the
transpose (D_a -> D_a^T on every copy), so that the stored integer tensors are exactly invariant.
Stored denominator: L' = 2^log2L * 2 * n!.

usage: python make_cert_sym.py sol.npz n a_num a_den t_num t_den log2L out.txt
"""
import sys
import math
import itertools
import numpy as np
import brauer as br

sol, n = sys.argv[1], int(sys.argv[2])
anum, aden, tnum, tden = map(int, sys.argv[3:7])
log2L = int(sys.argv[7])
out = sys.argv[8]
d = np.load(sol)
L = 2 ** log2L
tr = br.pt_perm([0, 1, 2])          # full transpose of a diagram


def exact_sym(x):
    X = np.rint(x * L).astype(np.int64).astype(object).reshape([br.NDIAG] * n)
    acc = np.zeros_like(X)
    for perm in itertools.permutations(range(n)):
        Y = np.transpose(X, perm)
        acc = acc + Y
        # transpose on every copy: Y'[tr(a1),...,tr(an)] = Y[a1,...,an]
        Yt = np.empty_like(Y)
        for idx in itertools.product(range(br.NDIAG), repeat=n):
            Yt[tuple(tr[i] for i in idx)] = Y[idx]
        acc = acc + Yt
    return acc.reshape(-1)


XA = exact_sym(d['xA'])
XB = exact_sym(d['xB'])
Lp = L * 2 * math.factorial(n)
print('max |xA|', np.abs(d['xA']).max(), 'max |xB|', np.abs(d['xB']).max())
print('symmetrisation change (max, in units of 1/L):',
      np.abs(np.array(XA, dtype=float) / (2 * math.factorial(n)) - d['xA'] * L).max(),
      np.abs(np.array(XB, dtype=float) / (2 * math.factorial(n)) - d['xB'] * L).max())
with open(out, 'w') as f:
    f.write(f'n {n}\n')
    f.write(f'a {anum} {aden}\n')
    f.write(f't {tnum} {tden}\n')
    f.write(f'L {Lp}\n')
    for k, m in enumerate(br.DIAGRAMS):
        f.write(f'D{k} ' + ' '.join(f'{p} {q}' for p, q in m) + '\n')
    f.write('QA ' + ' '.join(str(int(v)) for v in XA) + '\n')
    f.write('QB ' + ' '.join(str(int(v)) for v in XB) + '\n')
print('wrote', out)
