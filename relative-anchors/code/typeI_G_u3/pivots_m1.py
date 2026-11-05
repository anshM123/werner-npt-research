"""u = 3, m = 1: the LDL^T pivots of G_z (closed_m1_u3.pkl) as reduced rational functions over Q(s0, ..., s3).

G_z > 0 iff its four pivots p_1 = G_00, p_2, p_3, p_4 (symmetric Gaussian elimination in the order rho, z_1, z_2, z_3)
are positive.  Each pivot is computed in the rational function field (python-flint, gcd-reduced), and its numerator
and denominator are factored.  Output: pivots_m1.pkl = list of (num dict, den dict) over (s0, s1, s2, s3)."""
import sys, os, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from flintfield import Field, F as FF

u = 3
t0 = time.time()
D = pickle.load(open(f'closed_m1_u{u}.pkl', 'rb'))
K = Field(D['names'])
G = [[FF(K, K.ctx.from_dict(n_), K.ctx.from_dict(d_)) for (n_, d_) in row] for row in D['G']]
d = len(G)
A = [row[:] for row in G]
piv = []
for p in range(d):
    pv = A[p][p]
    piv.append(pv)
    for i in range(p + 1, d):
        f = A[i][p] / pv
        for j in range(p + 1, d):
            A[i][j] = A[i][j] - f * A[p][j]
    print(f'pivot {p + 1}: num deg {pv.n.total_degree()} ({len(pv.n)} terms), den deg {pv.d.total_degree()} ({len(pv.d)} terms) ({time.time() - t0:.1f} s)', flush=True)
out = []
for k, pv in enumerate(piv):
    cn, fn = pv.n.factor()
    cd, fd = pv.d.factor()
    print(f'pivot {k + 1}: num = {cn} * ' + ' * '.join(f'[deg {f.total_degree()}, {len(f)} t]^{e}' for f, e in fn), flush=True)
    print(f'          den = {cd} * ' + ' * '.join(f'[deg {f.total_degree()}, {len(f)} t]^{e}' for f, e in fd), flush=True)
    out.append((pv.n.to_dict(), pv.d.to_dict()))
pickle.dump(dict(pivots=out, names=D['names']), open('pivots_m1.pkl', 'wb'))
print(f'total {time.time() - t0:.1f} s')
print('PASS')
