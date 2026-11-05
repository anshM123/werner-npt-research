"""The two implementations of the equal-U-weight reduction agree: equ_sym.py (written for one level only; its output
logs/equ_sym.pkl) and lev.py with l = 1, offset 2 give the same rational functions F_inv (3 x 3) and F_std (2 x 2) in
(p, q, x), compared as elements of Q(p, q, x) (identical numerators and denominators after reduction).
Usage: python cmp_impl.py (after equ_sym.py)"""
import pickle
from flintfield import Field, F
from lev import build

D = pickle.load(open('logs/equ_sym.pkl', 'rb'))
S = build(1, [2], 1, want_phi=True)
K = S['K']                     # names p, q1, x1: same order as p, q, x
ctx = K.ctx


def rf(nd):
    n, d = nd
    return F(K, ctx.from_dict(n), ctx.from_dict(d))


ok = True
for i in range(3):
    for j in range(3):
        ok &= rf(D['Finv'][i][j]) == S['Finv'][i][j]
for i in range(2):
    for j in range(2):
        ok &= rf(D['Fstd'][i][j]) == S['Fstd'][0][i][j]
print(f'F_inv (9 entries) and F_std (4 entries) identical in both implementations: {ok}')
print('PASS' if ok else 'FAIL')
