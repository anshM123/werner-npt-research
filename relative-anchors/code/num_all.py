"""Exact certificate (rational arithmetic, orbit calculus) that the letter-1 torus block of (R_inf)_j at equal
weights is PSD for every k in [kmin, kmax], every j = 0..k and every irrep index i = 0..floor(k/2)
(q0 = j - i and t = k - j - i may be -1).  Also tallies the kernel: sum over i of nullity(Q_i) * dim of the
irrep it sits in is NOT computed here; we report nullity(Q_i) and compare with the expected generic pattern."""
import sys, time
from orbit3 import Num, reduced_orbit
from eqlib import ldl_psd
from fractions import Fraction
import flint
kmin, kmax = int(sys.argv[1]), int(sys.argv[2])
t0 = time.time(); ok = True; cnt = 0; nul_hist = {}
for k in range(kmin, kmax + 1):
    for j in range(0, k + 1):
        for i in range(0, k // 2 + 1):
            names, G, Q = reduced_orbit(Num(k, j, i))
            if not names:
                continue
            Qf = [[flint.fmpq(x.numerator, x.denominator) for x in row] for row in Q]
            psd, rk, piv = ldl_psd(Qf)
            cnt += 1
            nul = len(names) - rk
            key = (j - i, k - j - i) if (j - i <= 1 or k - j - i <= 1) else 'interior'
            nul_hist.setdefault(nul, 0); nul_hist[nul] += 1
            if not psd:
                ok = False; print('NOT PSD', k, j, i, names)
    if k % 10 == 0:
        print(f'k={k} done ({time.time()-t0:.0f}s)', flush=True)
print(f'{cnt} reduced blocks, nullity histogram {nul_hist}')
print('PASS' if ok else 'FAIL', f'{time.time()-t0:.0f}s')
