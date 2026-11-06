"""The third LDL^T pivot of (B_inv) for u1, u2 >= 2 (s_t = s or 1), computed after composing the 3 x 3 diagonal block
into the orthant variables (x1, x2, r, s) (u_i = 2 + x_i, sigma = s_t + r, s1 = s); same certificate as prove_gersh.py.
Usage: python b3early.py s|1 [bz]"""
import sys, os, time, gc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
import prove_gersh as PG
from flintfield import F as FE

Q = flint.fmpq
ka = sys.argv[1]
mode = 'T11' if ka == 's' else 'T22'
BZ = 'bz' in sys.argv[2:]
KX, GX, ONE = PG.KX, PG.GX, PG.ONE
x1, x2, r, s = GX['x1'], GX['x2'], GX['r'], GX['s']
F = PG.F
t0 = time.time()
names, keep, G = PG.load('tl_inv_gen_%s.pkl' % mode)
a = F(s) if ka == 's' else F(ONE)
sig = a + F(r)
mp = dict(u1=2 * ONE + x1, u2=2 * ONE + x2, sig=sig.n, s1=s)
Gz = [[PG.comp(G[(min(i, j), max(i, j))], names, mp) for j in range(3)] for i in range(3)]
print('composed: entry sizes %s [%.1fs]' % ([len(Gz[i][j].n) for i in range(3) for j in range(i, 3)], time.time() - t0), flush=True)
del G
gc.collect()
S1 = (2 * ONE + x1) * s
S2 = 2 * ONE + x2
S = S1 + S2
H = PG.H
one = F(ONE)
om_th = one * H * F(r) / sig
if BZ:
    c1 = a * sig / (F(S) * (sig + F(S)))
    c2 = a * F(S1 * S2) / F(S)
    W = [[one, F(0 * ONE), F(0 * ONE)], [F(0 * ONE), c1 * F(S1 * S1) + c2, c1 * F(S1 * S2) - c2],
         [F(0 * ONE), c1 * F(S1 * S2) - c2, c1 * F(S2 * S2) + c2]]
    A = [[Gz[i][j] - om_th * W[i][j] for j in range(3)] for i in range(3)]
else:
    Nm = [[one, F(0 * ONE), F(0 * ONE)], [F(0 * ONE), FE(KX, ONE, S), FE(KX, S2, S)], [F(0 * ONE), FE(KX, ONE, S), FE(KX, -S1, S)]]
    T = [[sum((Gz[k][l] * Nm[l][j] for l in range(3) if Nm[l][j] != 0), F(0 * ONE)) for j in range(3)] for k in range(3)]
    Gy = [[sum((Nm[k][i] * T[k][j] for k in range(3) if Nm[k][i] != 0), F(0 * ONE)) for j in range(3)] for i in range(3)]
    kap = [one, a * sig / (F(S) * (sig + F(S))), a * F(S1 * S2) / F(S)]
    A = [[Gy[i][j] - (om_th * kap[i] if i == j else F(0 * ONE)) for j in range(3)] for i in range(3)]
    del T, Gy
del Gz
gc.collect()
print('matrix built [%.1fs]' % (time.time() - t0), flush=True)
ok = True
for p in range(3):
    pv = A[p][p]
    ok &= PG.cert('tl_inv_gen B_inv s_t = %s pivot %d (orthant variables%s)' % (ka, p + 1, ', z-coordinates' if BZ else ''),
                  pv, ('x1', 'x2', 'r'))
    for i in range(p + 1, 3):
        f = A[i][p] / pv
        for j in range(p + 1, 3):
            A[i][j] = A[i][j] - f * A[p][j]
        gc.collect()
print('done [%.1fs]' % (time.time() - t0))
print('PASS' if ok else 'FAIL')
