"""Coverage of the low irreps for every k: for every k = 2..KCHK, every j = -1..k and every i <= B with a nonempty
twisted space V_i, either k <= KFIX (fixed-k certificates) or some proved family contains (q0, t, i) (a fixed entry
must match; a symbolic entry needs q0 >= 4, t >= 4, i >= 2).
Usage: python coverage_irr.py KCHK KFIX B results1.pkl [results2.pkl ...]"""
import sys, pickle, itertools
from fractions import Fraction as Fr
import mlib
from famlib import SMIN

KCHK, KFIX, B = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
res = {}
for f in sys.argv[4:]:
    res.update(pickle.load(open(f, 'rb')))
proved = {key for key, r in res.items() if r.startswith('PD')}
ok = True
nfam = nfix = 0
for k in range(2, KCHK + 1):
    u = k - 1
    for j in range(-1, k + 1):
        for i in range(0, min(B, u // 2) + 1):
            C = mlib.Calc(k, j, i, sigma=Fr(1))
            if not any(mlib.norm_of(C, nm) for nm in mlib.NAMES):
                continue
            if k <= KFIX:
                nfix += 1
                continue
            q0, t = j - i, u - j - i
            cands = itertools.product((q0, None) if q0 >= SMIN['q0'] else (q0,),
                                      (t, None) if t >= SMIN['t'] else (t,),
                                      (i, None) if i >= SMIN['i'] else (i,))
            if any(c in proved for c in cands):
                nfam += 1
            else:
                ok = False
                print('NOT COVERED', k, j, i, (q0, t, i))
print(f'k = 2..{KCHK}, every j, i <= {B}: {nfam} nonempty reduced forms with k > {KFIX} in proved families, '
      f'{nfix} with k <= {KFIX}')
print('PASS' if ok else 'FAIL')
