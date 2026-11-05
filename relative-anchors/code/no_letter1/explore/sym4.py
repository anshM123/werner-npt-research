"""H for a two-atom measure m1 delta_{xT + s1} + m2 delta_{xT + s2}, e = m1 + m2: coefficient signs."""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import flint
from flintfield import Field
t0 = time.time()
K = Field(['p', 'q', 'r', 'm1', 'm2', 's1', 's2'])
p, q, r, m1, m2, s1, s2 = K.gens
e = m1 + m2
xT = p + q + r; c = xT + 2 * e
gi = p / (c - p); gl = q / (c - q); Gi = p / (c - p + q); Gl = q / (c - q + p)
dl = (1 + gi) * (1 + gl); ep = (1 + Gi) * (1 + Gl)
t1 = xT + s1; t2 = xT + s2
Ph = lambda d: m1 / (t1 + d) + m2 / (t2 + d)
X = 2 + Ph(0); Y = 2 + Gi + Gl + Ph(p + q); Z = 2 + gi + Ph(p); W = 2 + gl + Ph(q)
D = dl * X * Y - ep * Z * W
print('D built', len(D.n), f'[{time.time()-t0:.1f}s]', flush=True)
negD = [(ex, cf) for ex, cf in D.n.to_dict().items() if cf < 0]
print('D two atoms: num terms', len(D.n), 'negative', len(negD), flush=True)
Fii = 2 - X * (1 + gi) * (1 + gi) / (Z * Z)
H = q * Fii - (q + r) * D / (Z * W * Y)
n = H.n
neg = [(ex, cf) for ex, cf in n.to_dict().items() if cf < 0]
cst, fac = H.d.factor()
print('H two atoms: num terms', len(n), 'deg', n.total_degree(), 'negative', len(neg), 'den factors positive', all(all(cf > 0 for cf in f.coeffs()) for f, _ in fac), f'[{time.time()-t0:.1f}s]')
print(neg[:10])
