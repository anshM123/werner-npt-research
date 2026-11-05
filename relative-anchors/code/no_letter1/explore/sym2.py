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
    neg = [(ex, cf) for ex, cf in n.to_dict().items() if cf < 0]
    cst, fac = d.factor()
    dpos = all(all(cf > 0 for cf in f.coeffs()) for f, _ in fac) and cst > 0
    print(f'{tag}: num {len(n)} terms deg {n.total_degree()}, negative coeffs {len(neg)}; den positive factors {dpos} [{time.time()-t0:.1f}s]', flush=True)
    return neg
def Kfull(t, s):
    K0 = dl * (2 + Gi + Gl) / (t * t) + 2 * dl / ((t + a + b) * (t + a + b)) - ep * (2 + gl) / ((t + a) * (t + a)) - ep * (2 + gi) / ((t + b) * (t + b))
    K1 = dl / ((s + a + b) * t * t) + dl / (s * (t + a + b) * (t + a + b)) - ep / ((s + b) * (t + a) * (t + a)) - ep / ((s + a) * (t + b) * (t + b))
    return K0, K1
K0, K1 = Kfull(xT + tau, xT + sig)
neg = info(K0 + e * K1, 'K0 + e K1 shifted')
print(neg[:20])
# s in [xT, c]: s = xT + 2 e sig/(1 + sig)  (sig >= 0 maps onto [xT, c))
for nm, t, s in (('t shifted, s compact', xT + tau, xT + 2 * e * sig / (1 + sig)),
                 ('both compact', xT + 2 * e * tau / (1 + tau), xT + 2 * e * sig / (1 + sig))):
    K0, K1 = Kfull(t, s)
    info(K1, 'K1 ' + nm)
    info(K0 + e * K1, 'K0 + e K1 ' + nm)
