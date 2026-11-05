import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import flint
from flintfield import Field
t0 = time.time()
K = Field(['p', 'q', 'r', 'e', 'lam'])
p, q, r, e, lam = K.gens
xT = p + q + r; c = xT + 2 * e
gi = p / (c - p); gl = q / (c - q); Gi = p / (c - p + q); Gl = q / (c - q + p)
dl = (1 + gi) * (1 + gl); ep = (1 + Gi) * (1 + Gl)
def info(R, tag):
    n = R.n; d = R.d
    neg = [(ex, cf) for ex, cf in n.to_dict().items() if cf < 0]
    cst, fac = d.factor()
    dpos = all(all(cf > 0 for cf in f.coeffs()) for f, _ in fac) and cst > 0
    cn, fn = n.factor()
    print(f'{tag}: num {len(n)} terms deg {n.total_degree()}, negative coeffs {len(neg)}; num factors {[(len(f), f.total_degree(), m) for f, m in fn]}; den positive factors {dpos} [{time.time()-t0:.1f}s]', flush=True)
    return neg
def Hof(tau):
    Ph = lambda d: e / (tau + d)
    X = 2 + Ph(0); Y = 2 + Gi + Gl + Ph(p + q); Z = 2 + gi + Ph(p); W = 2 + gl + Ph(q)
    D = dl * X * Y - ep * Z * W
    Fii = 2 - X * (1 + gi) * (1 + gi) / (Z * Z)
    return q * Fii - (q + r) * D / (Z * W * Y), Fii, D
H1, Fii1, D1 = Hof(xT)
neg = info(H1, 'H one atom (|U| = 1)'); print(neg[:10])
info(Fii1, 'Fii one atom')
Hc, Fiic, Dc = Hof(c)
neg = info(Hc, 'H tiny limit (e delta_c)'); print(neg[:10])
Hl, _, _ = Hof(xT + 2 * e * lam)
neg = info(Hl, 'H single atom at xT + 2 e lam'); print(len(neg))
