"""Test 6 (exact): blocks without the obstruction of newton_obstr.py but with a non-PSD coefficient matrix at q = 1.
Try the multipliers q = |s|^N, N = 0..NMAX (Polya-type), and q = products of the irreducible factors of den with
exponents 0..1: if every monomial coefficient matrix of q * den * M_g is PSD, M_g = sum_m (s^m / (q den)) C_m with C_m
PSD is a fixed-vector representation (LDL^T of each C_m).
Usage: python polya_mult.py k NMAX blockname [blockname ...]"""
import sys, itertools, time
from psdtools import load_blocks, poly_matrix, coeff_mats, psd_exact

k = int(sys.argv[1])
NMAX = int(sys.argv[2])
names = sys.argv[3:]
t0 = time.time()
for blk in load_blocks(k):
    if blk['name'] not in names:
        continue
    ctx = blk['ctx']
    one = ctx.from_dict({(0,) * len(blk['names']): 1})
    S = sum(ctx.gens(), ctx.from_dict({}))
    if k == 4:
        S = S + one
    den, Nt = poly_matrix(blk)
    cst, fac = den.factor()
    found = None
    q = one
    for N in range(NMAX + 1):
        if N > 0:
            q = q * S
        _, Ntq = poly_matrix(blk, q)
        C = coeff_mats(Ntq)
        nf = sum(1 for m in C if not psd_exact(C[m])[0])
        print(f'{blk["name"]}: q = |s|^{N}: {nf} non-PSD of {len(C)} coefficient matrices  [{time.time() - t0:.1f}s]',
              flush=True)
        if nf == 0:
            found = f'|s|^{N}'
            break
    if found is None:
        for es in itertools.product((0, 1), repeat=len(fac)):
            if not any(es):
                continue
            q = one
            for (f, e), x in zip(fac, es):
                if x:
                    q = q * f
            _, Ntq = poly_matrix(blk, q)
            C = coeff_mats(Ntq)
            nf = sum(1 for m in C if not psd_exact(C[m])[0])
            if nf == 0:
                found = f'product of den factors {es}'
                break
        print(f'{blk["name"]}: products of den factors: {"found " + found if found else "none works"}', flush=True)
    print(f'RESULT {blk["name"]}: {"representation found with q = " + found if found else "no representation found"}')
print('PASS')
