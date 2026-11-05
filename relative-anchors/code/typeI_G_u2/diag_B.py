import pickle, itertools
from fractions import Fraction as Fr
import numpy as np
res = pickle.load(open('build_B.pkl', 'rb'))
H = Fr(2, 5)
def ev(dct, x):
    return sum(float(c) * x[0] ** int(m[0]) * x[1] ** int(m[1]) * x[2] ** int(m[2]) for m, c in dct.items())
for name, dct, degs, rem in res[H]:
    # vertices
    vals = {v: ev(dct, v) for v in itertools.product((0., 1.), repeat=3)}
    print(name, 'vertices:', {k: f'{v:.3g}' for k, v in vals.items()})
    g = np.linspace(0, 1, 11)
    mn = min((ev(dct, (a, b, c)), (a, b, c)) for a in g for b in g for c in g)
    print('   grid min', mn)
    # edges th=1, mu=0 ; th=1, mu=1 ; etc
    for fixed in [(1, 0, None), (1, 1, None), (None, 0, 1), (None, 1, 1), (1, None, 1), (0, None, 0), (None, 0, 0), (1, None, 0)]:
        vs = []
        for t in g:
            x = [t if f is None else f for f in fixed]
            vs.append(ev(dct, x))
        print('   edge', fixed, 'min', f'{min(vs):.3g}', 'max', f'{max(vs):.3g}')
