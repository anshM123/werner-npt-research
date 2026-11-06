"""Coverage bookkeeping for every (k, j, i).  For every k = 2..KCHK, every j = -1..k and every i with a nonempty
twisted space V_i: covered if k <= KFIX (fixed-k certificates, Theorem M) or if some proved family contains
(q0, t, i) = (j - i, k - 1 - j - i, i) (a fixed entry must match; a symbolic entry needs q0 >= 4, t >= 4, i >= 2).
Also: every nonempty form has q0 >= -1 and t >= -2, and no family reported empty contains a nonempty form.
The claim checked: every nonempty form that is not covered lies in the open region REGION, either 'corner:B'
(q0 >= 4, t >= 4, i >= B) or 'band:B' (9 <= j <= k - 9, i >= B); the script also prints the smallest k of an uncovered
form and counts the covered forms of the middle region (k > KFIX, 9 <= j <= k - 9) by family.
A pickle named members*.pkl holds individually certified forms {(k, j, i): status} (cert_rest.py); those with status
'PD on complement' count as covered.
Usage: python coverage_mid.py KCHK KFIX REGION results1.pkl [results2.pkl ...]"""
import sys, pickle, itertools
from collections import Counter
from fractions import Fraction as Fr
import mlib
from famlib import SMIN

KCHK, KFIX = int(sys.argv[1]), int(sys.argv[2])
RKIND, BREST = sys.argv[3].split(':')
BREST = int(BREST)
kmin_unc = None
res = {}
members = set()
for f in sys.argv[4:]:
    d = pickle.load(open(f, 'rb'))
    if f.split('/')[-1].startswith('members'):
        members |= {key for key, r in d.items() if r.startswith('PD')}
    else:
        res.update(d)
proved = {key for key, r in res.items() if r.startswith('PD')}
empty = {key for key, r in res.items() if r == 'empty'}
ok = True
cnt = Counter()
bym = Counter()
uncovered = 0
for k in range(2, KCHK + 1):
    u = k - 1
    for j in range(-1, k + 1):
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
            if k <= KFIX:
                cnt['fixed'] += 1
                continue
            if (k, j, i) in members:
                cnt['member'] += 1
                continue
            fams = [c for c in cands if c in proved]
            if fams:
                cnt['family'] += 1
                if 9 <= j <= k - 9:
                    bym[fams[0]] += 1
                continue
            uncovered += 1
            kmin_unc = k if kmin_unc is None else min(kmin_unc, k)
            inside = (q0 >= 4 and t >= 4) if RKIND == 'corner' else (9 <= j <= k - 9)
            if not (inside and i >= BREST):
                ok = False
                print('NOT COVERED AND NOT IN THE OPEN CORNER', k, j, i, (q0, t, i))
print(f'k = 2..{KCHK}, every j and i: {cnt["fixed"]} nonempty forms with k <= {KFIX}, {cnt["family"]} in proved '
      f'families, {cnt["member"]} certified one by one (k > {KFIX}), {uncovered} not covered, all in the region '
      f'{sys.argv[3]} (smallest k of an uncovered form: {kmin_unc})')
print('middle region (k > %d, 9 <= j <= k - 9), forms per family:' % KFIX)
lab = lambda v, nm: f'{nm}>={SMIN[nm]}' if v is None else f'{nm}={v}'
for fam, c in sorted(bym.items(), key=lambda z: str(z[0])):
    print(f'  [{lab(fam[0], "q0")}, {lab(fam[1], "t")}, {lab(fam[2], "i")}]: {c}')
print('PASS' if ok else 'FAIL')
