"""Coverage by the family certificates.  For every k = 2..KCHK and every j in the claimed region (REGION = 'band:J:D':
j <= J or k - j <= D), every i with a nonempty twisted space V_i must be covered: k <= KFIX (the fixed-k certificates),
or some proved family contains (q0, t, i) = (j - i, k - 1 - j - i, i) (a fixed entry must match; a symbolic entry needs
the value >= its minimum: q0 >= 4, t >= 4, i >= 2).  Also: every nonempty (k, j, i) has q0 >= -1 and t >= -2, and no
family reported 'empty' contains a nonempty (k, j, i) with k <= KCHK.
Usage: python coverage_fam.py KCHK KFIX REGION results1.pkl [results2.pkl ...]"""
import sys, pickle, itertools
from fractions import Fraction as Fr
import mlib
from famlib import SMIN

KCHK, KFIX, REGION = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
res = {}
for f in sys.argv[4:]:
    res.update(pickle.load(open(f, 'rb')))
_, J, D = REGION.split(':')
J, D = int(J), int(D)
proved = {key for key, r in res.items() if r.startswith('PD')}
empty = {key for key, r in res.items() if r == 'empty'}
ok = True
cnt = {'family': 0, 'fixed': 0}
for k in range(2, KCHK + 1):
    u = k - 1
    for j in range(-1, k + 1):
        inregion = j <= J or k - j <= D
        for i in range(0, u // 2 + 1):
            C = mlib.Calc(k, j, i, sigma=Fr(1))
            nonempty = any(mlib.norm_of(C, nm) for nm in mlib.NAMES)
            q0, t = j - i, u - j - i
            cands = list(itertools.product((q0, None) if q0 >= SMIN['q0'] else (q0,),
                                           (t, None) if t >= SMIN['t'] else (t,),
                                           (i, None) if i >= SMIN['i'] else (i,)))
            if nonempty and any(c in empty for c in cands):
                ok = False
                print('FAMILY REPORTED EMPTY BUT CONTAINS A NONEMPTY FORM', k, j, i)
            if not nonempty:
                continue
            if q0 < -1 or t < -2:
                ok = False
                print('NONEMPTY OUTSIDE THE PARAMETER RANGE', k, j, i)
            if not inregion:
                continue
            if k <= KFIX:
                cnt['fixed'] += 1
                continue
            if any(c in proved for c in cands):
                cnt['family'] += 1
            else:
                ok = False
                print('NOT COVERED', k, j, i, (q0, t, i))
print(f"k = 2..{KCHK}, j <= {J} or k - j <= {D}: {cnt['family']} nonempty reduced forms with k > {KFIX} in proved "
      f"families, {cnt['fixed']} with k <= {KFIX} (fixed-k certificates); {len(proved)} proved families")
print('PASS' if ok else 'FAIL')
