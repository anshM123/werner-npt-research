"""Negative controls (exact): the vector v of Conjecture M cannot be changed much.  With beta doubled
(v' = (4 s_t (s_t + S), 6 (2 s_t + S), ...)) or halved (v'' = (4 s_t (s_t + S), (3/2)(2 s_t + S), ...)), G_z v' (resp.
G_z v'') has a negative entry at explicit rational points, for u = 1, 3, 5 (equal and unequal U-weights).  The original
v passes at the same points.  Also, at equal U-weights the symbolic (D) certificate of check_equal_u.py fails (mixed
coefficients) for these perturbed vectors, as it must."""
import sys, os, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
import sympy as sp
from sympy import *  # noqa: F401,F403
from gzexact import Gz_exact

Q = flint.fmpq
ok = True


def Gv(G, v):
    n = G.nrows()
    return [sum((G[i, j] * v[j] for j in range(n)), Q(0)) for i in range(n)]


cases = [
    ('beta x 2', Q(2), Q(1), lambda u: [Q(1, 10)] * u),                       # S small: beta too large
    ('beta x 2', Q(2), Q(1), lambda u: [Q(1, 10 * (i + 1)) for i in range(u)]),
    ('beta x 1/2', Q(1, 2), Q(1), lambda u: [Q(10)] * u),                      # S large: beta too small
    ('beta x 1/2', Q(1, 2), Q(1), lambda u: [Q(10 + 3 * i) for i in range(u)]),
]
for name, fac, wt, wUf in cases:
    for u in (1, 3, 5):
        wU = wUf(u)
        G = Gz_exact(wt, wU)
        st = wt * wt
        S = sum((x * x for x in wU), Q(0))
        v = [4 * st * (st + S)] + [3 * (2 * st + S)] * u
        vp = [v[0]] + [fac * x for x in v[1:]]
        good_orig = all(x > 0 for x in Gv(G, v))
        bad_pert = any(x < 0 for x in Gv(G, vp))
        ok &= good_orig and bad_pert
        print(f'{name}, u = {u}, w_t = {wt}, w_U = {wU[:3]}{"..." if u > 3 else ""}: original v passes {good_orig}; '
              f'perturbed vector has a negative entry {bad_pert}', flush=True)

# symbolic: perturbed (D) at equal U-weights fails the coefficient test
u, sig, s = sp.symbols('u sigma s', positive=True)
x = sp.symbols('x', nonnegative=True)
R = {k: eval(v) for k, v in pickle.load(open('equal_u.pkl', 'rb')).items()}
Qinv = R['Qinv']
S = u * s
for fac, name in ((2, 'beta x 2'), (sp.Rational(1, 2), 'beta x 1/2')):
    v0 = 4 * sig * (sig + S)
    v1 = 3 * (2 * sig + S) * fac
    mixed = False
    for row in range(2):
        e = Qinv[row, 0] * v0 + Qinv[row, 1] * v1
        num, den = sp.fraction(sp.together(e))
        p = sp.Poly(sp.expand(num.subs(u, 1 + x)), sig, s, x)
        co = p.coeffs()
        if not (all(c > 0 for c in co) or all(c < 0 for c in co)):
            mixed = True
    ok &= mixed
    print(f'equal U-weights, {name}: the (D) numerators have mixed coefficients (certificate fails) {mixed}')
print('PASS' if ok else 'FAIL')
