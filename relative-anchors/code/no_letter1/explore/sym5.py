"""2-kernel of D (homogenised with the mass functional): coefficient signs with t = xT + s, t' = xT + s'."""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import flint
from flintfield import Field
t0 = time.time()
K = Field(['a', 'b', 'r', 'e', 's1', 's2'])
a, b, r, e, s1, s2 = K.gens
xT = a + b + r; c = xT + 2 * e
gi = a / (c - a); gl = b / (c - b); Gi = a / (c - a + b); Gl = b / (c - b + a)
dl = (1 + gi) * (1 + gl); ep = (1 + Gi) * (1 + Gl)
def psi(t):
    one = K(1)
    return (2 * one / e + one / t, (2 + Gi + Gl) / e + one / (t + a + b), (2 + gi) / e + one / (t + a), (2 + gl) / e + one / (t + b))
for nm, t, tp in (('free', xT + s1, xT + s2), ('compact', xT + 2 * e * s1 / (1 + s1), xT + 2 * e * s2 / (1 + s2))):
    X1, Y1, Z1, W1 = psi(t); X2, Y2, Z2, W2 = psi(tp)
    Dk = dl * (X1 * Y2 + X2 * Y1) - ep * (Z1 * W2 + Z2 * W1)
    n = Dk.n
    neg = [(ex, cf) for ex, cf in n.to_dict().items() if cf < 0]
    print(nm, 'D 2-kernel: num terms', len(n), 'negative', len(neg), f'[{time.time()-t0:.1f}s]')
    print(neg[:6])
d = n.to_dict()
print('coefficients of a b e^10 s1^i s2^j:')
for i in range(6):
    print([str(d.get((1, 1, 0, 10, i, j), 0)) for j in range(6)])
print('max degrees', n.degrees())
# group by (s1, s2) exponent: list all monomials with e-degree >= 9 and s-degrees (i,j) near (2,2)
for (i, j) in [(1, 2), (2, 1), (2, 3), (3, 2), (1, 3), (3, 1), (2, 2), (1, 1), (3, 3), (0, 4), (4, 0), (0, 2), (2, 0), (4, 4), (2, 4), (4, 2)]:
    terms = [(ex, str(cf)) for ex, cf in d.items() if ex[4] == i and ex[5] == j and ex[3] >= 9]
    print((i, j), terms[:8])
