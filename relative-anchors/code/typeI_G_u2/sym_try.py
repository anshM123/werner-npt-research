"""Timing test: symbolic Schur complement on the T-blocks, u = 2, m = 3, over QQ(w0, w1, w2, b) (w_U = (1, b))."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flintfield import Field
from cf import TypeI

K = Field(['w0', 'w1', 'w2', 'b'])
w0, w1, w2, b = K.gens
t0 = time.time()
ti = TypeI([w0, w1, w2, K(1), b], [0, 1, 2], F=K)
M = ti.build()
print('built', ti.n, f'{time.time() - t0:.1f}s', flush=True)
U = ti.U; r0 = U[0]
tb = []
for t in ti.T:
    tb.append(ti.idx[('R', t)])
    for h in U:
        tb.append(ti.idx[('Z', t, h)])
el = [i for i in range(ti.n) if i not in tb and not (ti.labels[i][0] == 'P' and ti.labels[i][1] == r0)]
order = el + tb
A = [[M[i][j] for j in order] for i in order]
ne = len(el)
for p in range(ne):
    t1 = time.time()
    pv = A[p][p]
    assert pv != 0
    for i in range(p + 1, len(order)):
        if A[i][p] == 0:
            continue
        f = A[i][p] / pv
        for j in range(p + 1, len(order)):
            if A[p][j] != 0:
                A[i][j] = A[i][j] - f * A[p][j]
    sz = max(len(A[i][j].n) + len(A[i][j].d) for i in range(p + 1, len(order)) for j in range(p + 1, len(order)))
    print('eliminated', ti.labels[order[p]], f'{time.time() - t1:.1f}s', 'max terms', sz, flush=True)
G = [[A[i][j] for j in range(ne, len(order))] for i in range(ne, len(order))]
import pickle
pickle.dump([[ (str(x.n), str(x.d)) for x in row] for row in G], open('sym_try_G.pkl', 'wb'))
print('total', f'{time.time() - t0:.1f}s')
