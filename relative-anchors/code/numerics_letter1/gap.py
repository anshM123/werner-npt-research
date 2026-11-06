"""Floating point: the smallest eigenvalue above the known kernel of the letter-1 block of (R_inf)_j, in units 2^Z Y, at
two weight values (n1 copies of squared weight 1, n2 of squared weight q).  python gap.py KMIN KMAX"""
import sys, time
import numpy as np
from l1float import L1Block
def gap(n1, n2, q, j):
    s = [1.0] * n1 + [q] * n2
    B = L1Block(s, j)
    Y = B.build()
    ev = np.linalg.eigvalsh(Y)
    nk = B.dim(j) + (1 if j == len(s) - 2 else 0)
    return ev[nk - 1], ev[nk], ev[0]
for k in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
    for j in sorted(set([k // 2, (k - 1) // 2 + 1, 3])):
        if not (3 <= j <= k - 2):
            continue
        for n2 in sorted(set([1, k // 2, k - 2])):
            n1 = k - n2
            row = []
            for q in [1.0, 0.5, 0.25, 0.1, 0.01]:
                lastk, g, m = gap(n1, n2, q, j)
                row.append(f'{g:.4f}')
                assert abs(lastk) < 1e-9 and m > -1e-9, (k, j, n1, n2, q, lastk, m)
            print(f'k={k} j={j} (n1,n2)=({n1},{n2}): gap*2^Z at q=1,1/2,1/4,1/10,1/100: ' + ' '.join(row), flush=True)
