"""The families [q0 >= 4, t >= 4, i = b] (two symbolic sizes and sigma) for the given values of b, proved as in
prove_fam.py.  With the families of prove_fam.py 1 1 and the fixed-k certificates they cover every (k, j) at the given
i: the isotypic parts of the S_{k-1}-irreps (k-1-i, i) and (k-2-i, i, 1).
Usage: python prove_lowirr.py b1 [b2 ...]"""
import sys, time, pickle
import prove_fam

bs = [int(x) for x in sys.argv[1:]]
T0 = time.time()
allok = True
results = {}
for b in bs:
    fam = (None, None, b)
    res, dt, info = prove_fam.prove(*fam, verbose=True)
    results[fam] = res
    allok &= res.startswith('PD') or res == 'empty'
    print(f'{fam}: {res}  ({dt:.1f}s)' + (f'  largest pivot {info["maxterms"]} terms' if info.get('maxterms') else ''),
          flush=True)
with open('logs/fam_lowirr_' + '_'.join(map(str, bs)) + '.pkl', 'wb') as fh:
    pickle.dump(results, fh)
print(f'{time.time() - T0:.0f}s')
print('PASS' if allok else 'FAIL')
