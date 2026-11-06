"""Exploration: factor structure of the closed forms of equal_m.py (no claim rests on it)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from emlib import load
R = load()
names = ('u', 'sigma', 's', 's1', 's2')
def show(nd, label):
    n, d = nd
    for part, nm in ((n, 'num'), (d, 'den')):
        c, fl = part.factor()
        desc = []
        for f, e in fl:
            degs = [max((m[i] for m in f.monoms()), default=0) for i in range(5)]
            neg = sum(1 for cc in f.coeffs() if cc < 0)
            desc.append('(%d terms, degs %s, neg %d)^%d' % (len(f), degs, neg, e))
        print(label, nm, c, ' * '.join(desc))
show(R['inv'][0][0], 'inv d rr')
show(R['inv'][0][1], 'inv d rz')
show(R['inv'][1][1], 'inv d zz')
show(R['inv'][0][2], 'inv o rr')
show(R['inv'][0][3], 'inv o rz')
show(R['inv'][1][2], 'inv o zr')
show(R['inv'][1][3], 'inv o zz')
show(R['std'][0][0], 'std d')
show(R['std'][0][1], 'std o')
