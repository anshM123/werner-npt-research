"""Summary of the certified families in the result pickles: per file the number of families by status, the number
of sanity points (member, sigma) compared with tlib.py, the largest entry and the largest pivot (terms of numerator
or denominator), the largest form, and for the determinant route the slice counts and the terms of N~.
python summary.py PKL [PKL ...]"""
import sys, pickle
from collections import Counter

ok = True
tot = Counter()
for p in sys.argv[1:]:
    d = pickle.load(open(p, 'rb'))
    res, info = d['results'], d['info']
    st = Counter('PD' if r.startswith('PD') else ('TIME' if r == 'TIME' else ('empty' if r == 'empty' else 'other'))
                 for r in res.values())
    pts = sum(inf.get('sanity_points', 0) for inf in info.values())
    ment = max([inf.get('form_terms', 0) or 0 for inf in info.values()] or [0])
    mpiv = max([inf.get('maxterms', 0) or 0 for inf in info.values()] or [0])
    nvec = max([len(inf.get('names', [])) for inf in info.values()] or [0])
    line = (f'{p}: {len(res)} families: ' + ', '.join(f'{a} {b}' for a, b in sorted(st.items())) +
            f'; sanity (member, sigma) points {pts}; largest form {nvec} vectors; largest entry {ment} terms; '
            f'largest pivot {mpiv} terms')
    if any('B' in inf for inf in info.values()):
        bs = sorted({inf['B'] for inf in info.values() if 'B' in inf})
        nt = max(inf.get('Nterms', 0) for inf in info.values())
        tm = sum(inf.get('time', 0) for inf in info.values())
        line += f'; determinant route: B in {bs} (B + 1 slices each), N~ up to {nt} terms, {tm:.0f}s in all'
    print(line)
    if st.get('other'):
        for key, r in res.items():
            if not (r.startswith('PD') or r in ('TIME', 'empty')):
                print('   not proved:', key, r)
    tot.update(st)
print('all files: ' + ', '.join(f'{a} {b}' for a, b in sorted(tot.items())))
print('PASS' if ok else 'FAIL')
