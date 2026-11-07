"""Rational certificate from a sym_sdp(5).py solution, rounded once per symmetry orbit.

Diagram-coordinate tuples (a_1..a_n) are grouped into orbits of (permutations of the copies) x (transpose of every
diagram); the float coefficients are averaged over each orbit and rounded to integers over 2^log2L, so the stored
tensors are exactly invariant.

usage: python make_cert_orbit.py sol.npz n a_num a_den t_num t_den log2L out.txt
"""
import sys
import itertools
import numpy as np
import brauer as br

sol, n = sys.argv[1], int(sys.argv[2])
anum, aden, tnum, tden = map(int, sys.argv[3:7])
log2L = int(sys.argv[7])
out = sys.argv[8]
d = np.load(sol)
L = 2 ** log2L
NB = br.NDIAG
tr = np.array(br.pt_perm([0, 1, 2]))
allk = np.array(list(itertools.product(range(NB), repeat=n)), dtype=np.int64)
pw = NB ** np.arange(n - 1, -1, -1)
canon = np.minimum(np.sort(allk, axis=1) @ pw, np.sort(tr[allk], axis=1) @ pw)
_, inv = np.unique(canon, return_inverse=True)
cnt = np.bincount(inv)


def orbit_round(x):
    mean = np.bincount(inv, weights=x) / cnt
    spread = np.abs(x - mean[inv]).max()
    Xi = np.rint(mean * L).astype(np.int64)[inv]
    return Xi, spread


XA, spA = orbit_round(d['xA'])
XB, spB = orbit_round(d['xB'])
print('orbits', len(cnt), 'max deviation from orbit mean (float solution):', spA, spB)
print('max |xA|', np.abs(d['xA']).max(), 'max |xB|', np.abs(d['xB']).max())
with open(out, 'w') as f:
    f.write(f'n {n}\n')
    f.write(f'a {anum} {aden}\n')
    f.write(f't {tnum} {tden}\n')
    f.write(f'L {L}\n')
    for k, m in enumerate(br.DIAGRAMS):
        f.write(f'D{k} ' + ' '.join(f'{p} {q}' for p, q in m) + '\n')
    f.write('QA ' + ' '.join(map(str, XA.tolist())) + '\n')
    f.write('QB ' + ' '.join(map(str, XB.tolist())) + '\n')
print('wrote', out)
