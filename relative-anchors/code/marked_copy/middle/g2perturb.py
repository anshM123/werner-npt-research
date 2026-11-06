"""A control for the vanishing test of stage A (not needed by the proof): copy the 103 slices of the column y = 0 to
OUT2, multiply the slice z = 51 by 1 + 1/1000, and run stage A there; the coefficients above z^60 must then be
nonzero (for the unchanged column they all vanish).  python g2perturb.py OUT OUT2"""
import sys, os, pickle, shutil
from fractions import Fraction as Fr
import g2interp

out, out2 = sys.argv[1], sys.argv[2]
os.makedirs(out2, exist_ok=True)
for z in range(103):
    shutil.copy(os.path.join(out, f's_0_{z}.pkl'), os.path.join(out2, f's_0_{z}.pkl'))
fn = os.path.join(out2, 's_0_51.pkl')
with open(fn, 'rb') as fh:
    d = pickle.load(fh)
d2 = {}
for e, (p, q) in d.items():
    v = Fr(p, q) * Fr(1001, 1000)
    d2[e] = (v.numerator, v.denominator)
with open(fn, 'wb') as fh:
    pickle.dump(d2, fh)
res = g2interp.stage_a(out2, 0, 102, 60)
with open(os.path.join(out2, 'col_0.pkl'), 'rb') as fh:
    c = pickle.load(fh)
print(f'perturbed column y = 0: nonzero coefficients above z^60: {c["bad"]} (of {42 * len(c["data"])})')
print('FAIL' if not res else 'PASS')
