"""Runs prove_split.prove_split on a list of families and saves the results (as prove_fam.py does) for the coverage
scripts.  A family is written q0,t,i with 'n' for a symbolic size (q0 >= 4, t >= 4, i >= 2), e.g. n,n,1.
Usage: python prove_mid.py OUT.pkl FAM [FAM ...] [--bscale=1001/1000]"""
import sys, time, pickle
from fractions import Fraction as Fr
from prove_split import prove_split, parse
from famlib import SMIN

out = sys.argv[1]
bscale = Fr(1)
fams = []
for a in sys.argv[2:]:
    if a.startswith('--bscale='):
        bscale = Fr(a.split('=')[1])
    else:
        fams.append(tuple(parse(s) for s in a.split(',')))
T0 = time.time()
allok = True
results = {}
nproved = nempty = 0
lab = lambda v, nm: f'{nm}>={SMIN[nm]}' if v is None else f'{nm}={v}'
for fam in fams:
    res, dt, info = prove_split(*fam, verbose=True, bscale=bscale)
    results[fam] = res
    ok = res.startswith('PD') or res == 'empty'
    allok &= ok
    nproved += res.startswith('PD')
    nempty += res == 'empty'
    print(f'[{lab(fam[0], "q0")}, {lab(fam[1], "t")}, {lab(fam[2], "i")}] {res}  ({dt:.1f}s)'
          + (f'  largest pivot {info["maxterms"]} terms' if info.get('maxterms') else ''), flush=True)
if bscale == 1:
    with open(out, 'wb') as fh:
        pickle.dump(results, fh)
print(f'{nproved} families proved, {nempty} empty, of {len(fams)}; {time.time() - T0:.0f}s')
print('PASS' if allok else 'FAIL')
