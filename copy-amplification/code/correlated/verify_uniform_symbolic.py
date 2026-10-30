"""Independent check of a dimension-uniform symbolic certificate file uniform_d_symbolic_weights_*.txt.

Written separately from uniform_d_cert.py (which produced the files).  For every listed generator
g = (G cores, Sig factors, O factors, X factors) with weight w_g(d), it checks, exactly in sympy:
  (1) sum_g w_g(d) dist_g(d) == target(d) identically in d, coordinate by coordinate, where
        target_k = C(n,k) p^k (1-p)^(n-k),  p = D_a (1+alpha) / (D_s (1-alpha) + D_a (1+alpha)),
        dist_g = convolution of the factor distributions
        G: (D_s, r0 D_a)/(D_s + r0 D_a),  O: (1, 0),  X = P^s + b P^a, b = (d+1)/(d-1): (1/2, 1/2),
        Sig = |Phi><Phi|_A (x) |Phi><Phi|_B on two sites: (D_s, 0, D_a)/d^2;
  (2) every w_g(d) >= 0 for all real d >= D0: numerator and denominator are factored, and each factor's real roots in
      [D0, oo) are counted exactly (sympy count_roots), with signs fixed at d = D0 + 1;
  (3) at most m cores per generator and n sites in total.
usage: verify_uniform_symbolic.py weights.txt n m r0_num r0_den alpha_num alpha_den D0
"""
import ast
import sys
from math import comb

import sympy as sp

path = sys.argv[1]
n, m = int(sys.argv[2]), int(sys.argv[3])
r0 = sp.Rational(int(sys.argv[4]), int(sys.argv[5]))
alpha = sp.Rational(int(sys.argv[6]), int(sys.argv[7]))
D0 = int(sys.argv[8])
d = sp.Symbol('d')

Ds, Da = d * (d + 1) / 2, d * (d - 1) / 2
FAC = {'G': [Ds / (Ds + r0 * Da), r0 * Da / (Ds + r0 * Da)], 'O': [sp.Integer(1), sp.Integer(0)],
       'X': [sp.Rational(1, 2), sp.Rational(1, 2)], 'Sig': [Ds / d ** 2, sp.Integer(0), Da / d ** 2]}
SIZE = {'G': 1, 'O': 1, 'X': 1, 'Sig': 2}

rows = []
with open(path) as fh:
    for line in fh:
        line = line.strip()
        if not line.startswith('{'):
            continue
        key, expr = line.split('}:', 1)
        g = ast.literal_eval(key + '}')
        rows.append((g, sp.sympify(expr.strip(), locals={'d': d})))
print(f'{path}: {len(rows)} generators')

total = [sp.Integer(0)] * (n + 1)
for g, w in rows:
    assert sum(SIZE[k] * v for k, v in g.items()) == n, g
    assert g.get('G', 0) <= m, g
    dist = [sp.Integer(1)]
    for name in ('G', 'Sig', 'O', 'X'):
        for _ in range(g.get(name, 0)):
            f = FAC[name]
            new = [sp.Integer(0)] * (len(dist) + len(f) - 1)
            for i, x in enumerate(dist):
                for j, y in enumerate(f):
                    new[i + j] += x * y
            dist = new
    assert len(dist) == n + 1
    for k in range(n + 1):
        total[k] += w * dist[k]

p = Da * (1 + alpha) / (Ds * (1 - alpha) + Da * (1 + alpha))
ok_identity = True
for k in range(n + 1):
    diff = sp.cancel(sp.together(total[k] - comb(n, k) * p ** k * (1 - p) ** (n - k)))
    if diff != 0:
        ok_identity = False
        print('  identity FAILS at k =', k)
print('  identity holds identically in d:', ok_identity)

ok_sign = True
for g, w in rows:
    num, den = sp.fraction(sp.together(w))
    for poly_expr in (num, den):
        P = sp.Poly(sp.expand(poly_expr), d)
        if P.degree() <= 0:
            if P.as_expr() < 0 and poly_expr is num:
                pass
            continue
        nroots = P.count_roots(D0, None)
        if nroots != 0:
            ok_sign = False
            print('  root in [D0, oo) for', g)
    val = w.subs(d, D0 + 1)
    if not val >= 0:
        ok_sign = False
        print('  negative weight at D0 + 1 for', g, float(val))
print(f'  every weight >= 0 on [{D0}, oo):', ok_sign)
print('PASS' if ok_identity and ok_sign else 'FAIL')
