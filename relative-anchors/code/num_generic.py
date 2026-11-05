"""Exact elimination (Num2) of the sym forms of the generic family q0 >= 2, t >= 4, i >= 2, tau >= 1 for all
k = q0 + t + 2i + tau in [kmin, kmax] (the only case I family whose symbolic elimination is not finished).
Usage: python num_generic.py kmin kmax [tmin=4] [taumin=1]  (the open corner is tmin = 6, taumin = 2)"""
import sys, time
import flint
from orbit2 import Num2, reduced2
from eqlib import ldl_psd
kmin, kmax = int(sys.argv[1]), int(sys.argv[2])
TMIN = int(sys.argv[3]) if len(sys.argv) > 3 else 4
TAUMIN = int(sys.argv[4]) if len(sys.argv) > 4 else 1
ok = True; cnt = 0; t0 = time.time()
for k in range(kmin, kmax + 1):
    for i in range(2, k):
        for tau in range(TAUMIN, k):
            for q0 in range(2, k):
                t = k - tau - 2 * i - q0
                if t < TMIN:
                    continue
                j = q0 + i + tau
                names, G, Q = reduced2(Num2(k, j, tau, i), 'sym')
                Qf = [[flint.fmpq(x.numerator, x.denominator) for x in row] for row in Q]
                psd, rk, piv = ldl_psd(Qf)
                cnt += 1
                if not psd or len(names) - rk != 1:
                    ok = False; print('FAIL', k, j, tau, i, len(names), rk, psd)
    print(f'k={k}: {cnt} generic forms so far ({time.time()-t0:.0f}s)', flush=True)
print(cnt, 'generic sym forms;', 'PASS (PSD, nullity exactly 1)' if ok else 'FAIL', f'{time.time()-t0:.0f}s')
