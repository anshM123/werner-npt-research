"""Families for the bands of j at both ends, symbolic in k.
For a fixed j only the irreps i <= j + 1 occur (q0 = j - i >= -1), and for a fixed k - j only i <= k - j + 1 occur
(t = k - 1 - j - i >= -2).  So the band j <= W needs only the families [q0 = a, t >= 4 symbolic, i = b] with
a + b <= W, and the band k - j <= W only the families [q0 >= 4 symbolic, t = c, i = b] with c + b <= W - 1; the members of
the band outside these families have k <= 2W + 5 and are covered by the fixed-k certificates.  Each family is proved as
in prove_fam.py (sanity at sample points, kernel identities, gauge, greedy elimination with one-signed pivots after
the shift, sigma free).
Usage: python prove_band.py bottom|top W"""
import sys, time, pickle
import prove_fam

mode, W = sys.argv[1], int(sys.argv[2])
fams = []
for b in range(0, W + 2):
    if mode == 'bottom':
        for a in range(-1, W - b + 1):
            fams.append((a, None, b))
    else:
        for c in range(-2, W - b):
            fams.append((None, c, b))
T0 = time.time()
allok = True
results = {}
nproved = nempty = 0
for fam in fams:
    res, dt, info = prove_fam.prove(*fam)
    results[fam] = res
    ok = res.startswith('PD') or res == 'empty'
    allok &= ok
    nproved += res.startswith('PD')
    nempty += res == 'empty'
    print(f'{fam}: {res}  ({dt:.1f}s)' + (f'  largest pivot {info["maxterms"]} terms' if info.get('maxterms') else ''),
          flush=True)
with open(f'logs/band_{mode}_{W}.pkl', 'wb') as fh:
    pickle.dump(results, fh)
print(f'{mode} band, width {W}: {nproved} families proved, {nempty} empty, of {len(fams)}; {time.time() - T0:.0f}s')
print('PASS' if allok else 'FAIL')
