"""Exploration (no claim rests on it): the cost of the parts of one slice, and of two cheaper variants (entries not
reduced; the determinant with all entries over the common factor of their denominators).  Every variant must give the
same reduced fraction.  python em3prof.py s|1 c"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
import em3lib as L
from flintfield import Field, F as FE

ka, c = sys.argv[1], int(sys.argv[2])
A, (names, Ao) = L.build_A(ka, orig=True)
OS = L.OrigSlicer(ka, names, Ao)
N = 2 + c
t = time.time()
sp = {(i, j): (OS.spec_u1(Ao[i][j].n, N), OS.spec_u1(Ao[i][j].d, N)) for i in range(3) for j in range(i, 3)}
t_spec = time.time() - t
t = time.time()
M = [[None] * 3 for _ in range(3)]
for i in range(3):
    for j in range(i, 3):
        M[i][j] = FE(OS.K, *sp[(i, j)])
        M[j][i] = M[i][j]
t_red = time.time() - t
t = time.time()
d = L.det3(M)
t_det = time.time() - t
t = time.time()
dd = FE(L.K3, d.n.compose(*OS.back, ctx=L.K3.ctx), d.d.compose(*OS.back, ctx=L.K3.ctx))
t_back = time.time() - t
print('x1 = %d: specialisation %.1fs, entry gcds %.1fs, determinant %.1fs, composition back %.1fs; num %d, den %d terms'
      % (c, t_spec, t_red, t_det, t_back, len(dd.n), len(dd.d)), flush=True)
# variant 1: entries not reduced
t = time.time()
M1 = [[None] * 3 for _ in range(3)]
for i in range(3):
    for j in range(i, 3):
        M1[i][j] = FE(OS.K, *sp[(i, j)], reduce=False)
        M1[j][i] = M1[i][j]
d1 = L.det3(M1)
d1 = FE(OS.K, d1.n, d1.d)
t1 = time.time() - t
print('variant without entry gcds: determinant %.1fs, equal %s' % (t1, d1.n == d.n and d1.d == d.d), flush=True)
# variant 2: common denominator factor of all entries taken out (det = det(C M)/C^3)
t = time.time()
C = sp[(0, 0)][1]
for k in sp:
    C = C.gcd(sp[k][1])
M2 = [[None] * 3 for _ in range(3)]
for i in range(3):
    for j in range(i, 3):
        n_, d_ = sp[(i, j)]
        M2[i][j] = FE(OS.K, n_, d_ / C, reduce=False)
        M2[j][i] = M2[i][j]
d2 = L.det3(M2)
d2 = d2 / FE(OS.K, C ** 3, OS.K.ctx.from_dict({(0, 0, 0): 1}), reduce=False)
t2 = time.time() - t
print('variant with the common factor C (%d terms) taken out: %.1fs, equal %s' % (len(C), t2, d2.n == d.n and d2.d == d.d),
      flush=True)
print('PASS')
