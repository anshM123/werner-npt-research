"""Exact elimination (rational arithmetic, orbit calculus Num2) of every CASE I letter-2 form (sym and anti, every
i, tau, j) for k in [kmin, kmax]; independent of eqlib2/run2 (different code path), cost independent of k.
Usage: python num_letter2.py kmin kmax"""
import sys, time
import flint
from orbit2 import Num2, reduced2
from eqlib import ldl_psd
kmin, kmax = int(sys.argv[1]), int(sys.argv[2])
ok = True; cnt = 0; hist = {}; t0 = time.time()
for k in range(kmin, kmax + 1):
    for j in range(0, k + 1):
        for tau in range(1, min(k, j + 1) + 1):
            for i in range(0, (k - tau) // 2 + 1):
                for mode in (('sym', 'anti') if tau >= 2 else ('sym',)):
                    names, G, Q = reduced2(Num2(k, j, tau, i), mode)
                    if not names:
                        continue
                    Qf = [[flint.fmpq(x.numerator, x.denominator) for x in row] for row in Q]
                    psd, rk, piv = ldl_psd(Qf)
                    cnt += 1
                    nul = len(names) - rk
                    hist[(mode, nul)] = hist.get((mode, nul), 0) + 1
                    if not psd:
                        ok = False; print('NOT PSD', k, j, tau, i, mode)
    print(f'k={k} done, {cnt} forms so far ({time.time()-t0:.0f}s)', flush=True)
print(cnt, 'case I forms; (mode, nullity) histogram', hist)
print('PASS' if ok else 'FAIL', f'{time.time()-t0:.0f}s')
