"""Exploration: structure of the symbolic G (u = 2, m = 3) in z-coordinates zeta_{t,h} = w_t w_h z_{t,h}."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from flintfield import Field, F as FF
from cf import TypeI

K = Field(['w0', 'w1', 'w2', 'b'])
w0, w1, w2, b = K.gens
wt = [w0, w1, w2]; wU = [K(1), b]
t0 = time.time()
ti = TypeI(wt + wU, [0, 1, 2], F=K)
M = ti.build()
U = ti.U; r0 = U[0]
tb = []; scale = []
for t in ti.T:
    tb.append(ti.idx[('R', t)]); scale.append(K(1))
    for j, h in enumerate(U):
        tb.append(ti.idx[('Z', t, h)]); scale.append(wt[t] * wU[j])
el = [i for i in range(ti.n) if i not in tb and not (ti.labels[i][0] == 'P' and ti.labels[i][1] == r0)]
order = el + tb
A = [[M[i][j] for j in order] for i in order]
ne = len(el)
for p in range(ne):
    pv = A[p][p]
    for i in range(p + 1, len(order)):
        if A[i][p] == 0:
            continue
        f = A[i][p] / pv
        for j in range(p + 1, len(order)):
            if A[p][j] != 0:
                A[i][j] = A[i][j] - f * A[p][j]
G = [[A[ne + i][ne + j] * scale[i] * scale[j] for j in range(9)] for i in range(9)]
print('done', f'{time.time() - t0:.1f}s', flush=True)
def even(p):
    return all(all(e % 2 == 0 for e in m) for m in p.to_dict().keys())
for (i, j) in [(0, 0), (0, 1), (1, 1), (1, 2), (0, 3), (0, 4), (1, 4), (2, 5), (1, 5)]:
    x = G[i][j]
    print((i, j), 'num terms', len(x.n), 'den terms', len(x.d), 'even num/den', even(x.n), even(x.d), flush=True)
# factor a denominator
x = G[0][3]
fac = x.d.factor()
print('den factors of G[0][3]:')
for f_, e_ in fac[1]:
    print('  ', e_, len(f_), str(f_)[:200])
x = G[0][0]
fac = x.d.factor()
print('den factors of G[0][0]:')
for f_, e_ in fac[1]:
    print('  ', e_, len(f_), str(f_)[:200])
