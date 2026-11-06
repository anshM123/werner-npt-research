"""Coverage bookkeeping for every (k, j, i) with k <= KCHK, with the families of prove_gen.py (symbolic sizes with
explicit lower bounds, and certificates of the whole form, of the (u - i, i) part only, or of the (u - i - 1, i, 1)
part only) besides the older results.

Inputs: pickles of three kinds.
  - members*.pkl: individually certified forms {(k, j, i): status} (cert_rest.py / cert_rest_mi.py); status
    'PD on complement' covers the whole form;
  - prove_gen.py outputs {'results': {key: status}, 'parts': '32' | '2' | '3'}: key (q0, t, i) with a fixed value or
    ('>=', minimum); status 'PD..' / 'PART2 PD..' / 'PART3 PD..' certifies the parts named, 'empty' says V_i = 0;
  - older family results {(q0, t, i): status} with None for a symbolic size (q0 >= 4, t >= 4, i >= 2): whole form.
A form is covered if k <= KFIX (Theorem M), or a members pickle certifies it, or a whole-form family contains it, or
a part-2 family and a part-3 family (or a whole-form one) contain it.  Its three-row part (the S vectors, present
when some block of Lemma S is split) is covered if the form is covered or a part-3 family contains it.
Checked: every nonempty form has q0 >= -1, t >= -2; no family reported empty contains a nonempty form; every
nonempty form that is not covered lies in the declared open region: k >= KOPEN, 9 <= j <= k - 9, q0 >= 4, t >= 4
and i >= IOPEN; and the three-row part of every nonempty form is covered.
Usage: python coverage_mi.py KCHK KFIX IOPEN KOPEN pickles..."""
import sys, pickle
from collections import Counter
from fractions import Fraction as Fr
import mlib

KCHK, KFIX, IOPEN, KOPEN = (int(x) for x in sys.argv[1:5])
OLDMIN = (4, 4, 2)
whole, part2, part3, empty = [], [], [], []
members = set()


def norm_key(key):
    out = []
    for v, mn in zip(key, OLDMIN):
        out.append(('>=', mn) if v is None else v)
    return tuple(out)


for f in sys.argv[5:]:
    d = pickle.load(open(f, 'rb'))
    base = f.replace('\\', '/').split('/')[-1]
    if base.startswith('members'):
        members |= {key for key, r in d.items() if r.startswith('PD')}
        continue
    if isinstance(d, dict) and 'results' in d:
        parts = d['parts']
        for key, r in d['results'].items():
            if r == 'empty':
                empty.append(key)
            elif r.startswith('PD') or r.startswith('PART'):
                if parts == '32':
                    whole.append(key)
                elif parts == '2':
                    part2.append(key)
                elif parts == '3':
                    part3.append(key)
        continue
    for key, r in d.items():
        if r == 'empty':
            empty.append(norm_key(key))
        elif r.startswith('PD'):
            whole.append(norm_key(key))


def contains(key, q0, t, i):
    for spec, v in zip(key, (q0, t, i)):
        if isinstance(spec, tuple):
            if v < spec[1]:
                return False
        elif spec != v:
            return False
    return True


MINS = [sorted({f[c][1] for f in whole + part2 + part3 + empty if isinstance(f[c], tuple)}) for c in range(3)]
WHOLE, PART2, PART3, EMPTY = set(whole), set(part2), set(part3), set(empty)


def hit(S, q0, t, i):
    """is (q0, t, i) in some family of S?  Every key containing it has each entry equal to the value or symbolic
    with a minimum <= the value; those keys are enumerated and looked up"""
    opts = [[v] + [('>=', mn) for mn in MINS[c] if v >= mn] for c, v in enumerate((q0, t, i))]
    return any((a, b, c) in S for a in opts[0] for b in opts[1] for c in opts[2])


ok = True
cnt = Counter()
uncovered = 0
kmin_unc = None
bykind = Counter()
for k in range(2, KCHK + 1):
    u = k - 1
    for j in range(-1, k + 1):
        for i in range(0, u // 2 + 1):
            C = mlib.Calc(k, j, i, sigma=Fr(1))
            names = [nm for nm in mlib.NAMES if mlib.norm_of(C, nm)]
            nonempty = bool(names)
            q0, t = j - i, u - j - i
            if nonempty and hit(EMPTY, q0, t, i):
                ok = False
                print('FAMILY REPORTED EMPTY BUT CONTAINS A NONEMPTY FORM', k, j, i)
            if not nonempty:
                continue
            if q0 < -1 or t < -2:
                ok = False
                print('NONEMPTY OUTSIDE THE PARAMETER RANGE', k, j, i)
            has3 = any(all(nm in names for nm in ('A%d%d' % (x, y0), 'B%d0%d' % (x, y0), 'B%d1%d' % (x, y0)))
                       for x in (0, 1) for y0 in (0, 1))
            w = hit(WHOLE, q0, t, i)
            p2 = hit(PART2, q0, t, i)
            p3 = hit(PART3, q0, t, i)
            if k % 10 == 0 and j == k // 2:
                # spot check of the indexed lookup against the direct test
                assert w == any(contains(f, q0, t, i) for f in whole)
                assert p2 == any(contains(f, q0, t, i) for f in part2)
                assert p3 == any(contains(f, q0, t, i) for f in part3)
            if k <= KFIX:
                cnt['fixed'] += 1
                covered = True
            elif (k, j, i) in members:
                cnt['member'] += 1
                covered = True
            elif w:
                cnt['family'] += 1
                covered = True
            elif p2 and (p3 or not has3):
                cnt['parts'] += 1
                covered = True
            else:
                covered = False
            if covered:
                continue
            if has3 and not p3:
                ok = False
                print('THREE-ROW PART NOT COVERED', k, j, i, (q0, t, i))
            uncovered += 1
            bykind[('three-row part covered by a family' if p3 else 'three-row part NOT covered') if has3 else 'no three-row part'] += 1
            kmin_unc = k if kmin_unc is None else min(kmin_unc, k)
            if not (k >= KOPEN and 9 <= j <= k - 9 and q0 >= 4 and t >= 4 and i >= IOPEN):
                ok = False
                print('NOT COVERED AND NOT IN THE OPEN REGION', k, j, i, (q0, t, i))
print(f'k = 2..{KCHK}, every j and i: {cnt["fixed"]} nonempty forms with k <= {KFIX}; {cnt["member"]} certified one '
      f'by one; {cnt["family"]} in whole-form families; {cnt["parts"]} by a (u - i, i) family together with a '
      f'(u - i - 1, i, 1) family; {uncovered} not covered (smallest k {kmin_unc}), all with k >= {KOPEN}, '
      f'9 <= j <= k - 9, q0 >= 4, t >= 4, i >= {IOPEN}; three-row parts: '
      f'({dict(bykind)})')
print(f'families read: {len(whole)} whole-form, {len(part2)} (u - i, i) part, {len(part3)} (u - i - 1, i, 1) part, '
      f'{len(empty)} empty; {len(members)} certified members')
print('PASS' if ok else 'FAIL')
