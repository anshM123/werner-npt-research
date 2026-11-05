import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import flint
from flintfield import Field
t0 = time.time()
K = Field(['a', 'b', 'r', 'e', 'tau', 'sig'])
a, b, r, e, tau, sig = K.gens
xT = a + b + r; S = xT + e; c = xT + 2 * e
gi = a / (c - a); gl = b / (c - b); Gi = a / (c - a + b); Gl = b / (c - b + a)
dl = (1 + gi) * (1 + gl); ep = (1 + Gi) * (1 + Gl)
def info(R, tag):
    n = R.n; d = R.d
    neg = [ (ex, cf) for ex, cf in n.to_dict().items() if cf < 0]
    cst, fac = d.factor()
    dpos = all(all(cf > 0 for cf in f.coeffs()) for f, _ in fac) and cst > 0
    cstn, facn = n.factor()
    print(f'{tag}: num {len(n)} terms deg {n.total_degree()}, negative coeffs {len(neg)}; num factors {[ (len(f), f.total_degree(), m) for f, m in facn]} const {cstn}; den positive factors {dpos} [{time.time()-t0:.1f}s]', flush=True)
    return neg
for name, t in (('t free', tau), ('t = xT + tau', xT + tau)):
    K0 = dl * (2 + Gi + Gl) / (t * t) + 2 * dl / ((t + a + b) * (t + a + b)) - ep * (2 + gl) / ((t + a) * (t + a)) - ep * (2 + gi) / ((t + b) * (t + b))
    info(K0, 'K0 ' + name)
for name, t, s in (('free', tau, sig), ('shifted', xT + tau, xT + sig)):
    K1 = dl / ((s + a + b) * (t * t)) + dl / (s * ((t + a + b) * (t + a + b))) - ep / ((s + b) * ((t + a) * (t + a))) - ep / ((s + a) * ((t + b) * (t + b)))
    neg = info(K1, 'K1 ' + name)
    if neg and len(neg) < 20: print(neg)
Dc = dl * (2 + e / c) * (2 + Gi + Gl + e / (c + a + b)) - ep * (2 + gi + e / (c + a)) * (2 + gl + e / (c + b))
info(Dc, 'D(e delta_c)')
D1 = dl * (2 + e / xT) * (2 + Gi + Gl + e / (xT + a + b)) - ep * (2 + gi + e / (xT + a)) * (2 + gl + e / (xT + b))
info(D1, 'D(e delta_xT)')
