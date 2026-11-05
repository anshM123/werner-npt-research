"""Try the multipliers (s0 + c s1)^N, c in {1,2,3,4,8,16}/{1,2,3,4,8,16}, N = 1..6, on the k = 2 block j = 1, label (1, 1, 0):
reports the smallest number of non-PSD coefficient matrices and prints the coefficient matrices at q = 1.
Usage: python linear_mult.py"""
import sys
sys.path.insert(0, '.')
from fractions import Fraction as Fr
import flint
from psdtools import load_blocks, poly_matrix, coeff_mats, psd_exact
blk = [b for b in load_blocks(2) if b['name'] == 'k2_j1_(1, 1, 0)'][0]
ctx = blk['ctx']
s0, s1 = ctx.gens()
best = None
for cn in [1, 2, 3, 4, 8, 16]:
    for cd in [1, 2, 3, 4, 8, 16]:
        c = flint.fmpq(cn, cd)
        for N in range(1, 7):
            q = (s0 + c * s1) ** N
            _, Nt = poly_matrix(blk, q)
            C = coeff_mats(Nt)
            nf = sum(1 for m in C if not psd_exact(C[m])[0])
            if best is None or nf < best[0]:
                best = (nf, cn, cd, N)
            if nf == 0:
                print('FOUND', cn, cd, N); print('PASS'); sys.exit()
print('best', best)
# show the den and the non-PSD monomials for q = 1
den, Nt = poly_matrix(blk)
C = coeff_mats(Nt)
print('den', den)
for m in sorted(C):
    ok, w = psd_exact(C[m])
    print(m, ok, [[str(x) for x in row] for row in C[m]] if not ok else '')
print('PASS')
