"""u = 3, m = 1: exact one-point check and negative controls for the closed form closed_m1_u3.pkl.

(P) G_z > 0 at (s0, s1, s2, s3) = (1/2, 1/4, 1/6, 1/12) (exact LDL^T pivots all positive).
(C1) The bound is sharp in the sense that no uniform spectral gap holds in the rescaled frame: with
     W = diag(1, s0 (diag(s_U) - s_U s_U^T/(s0 + S))) (the identity of the frame (rho, c z, zeta') in z-coordinates),
     G_z - h W is NOT positive definite for h = 1/10 at an explicit rational point with S small (exact witness).
(C2) Perturbing the closed form by -1/1000 in the (rho, rho) entry makes it indefinite at an explicit rational point
     (so positivity is not automatic from the shape of the entries); exact witness."""
import pickle, time, flint
from fractions import Fraction as Fr

u = 3
D = pickle.load(open(f'closed_m1_u{u}.pkl', 'rb'))
ctx = flint.fmpq_mpoly_ctx.get(tuple(D['names']), 'lex')
NUM = [[ctx.from_dict(n_) for (n_, d_) in row] for row in D['G']]
DEN = ctx.from_dict(D['G'][0][0][1])


def Gz(s):
    q = [flint.fmpq(x.numerator, x.denominator) for x in s]
    dv = DEN(*q)
    out = []
    for row in NUM:
        r = []
        for p in row:
            v = p(*q) / dv
            r.append(Fr(int(v.p), int(v.q)))
        out.append(r)
    return out


def ldl(A):
    A = [row[:] for row in A]
    n = len(A)
    piv = []
    for p in range(n):
        piv.append(A[p][p])
        if A[p][p] <= 0:
            return piv
        for i in range(p + 1, n):
            f = A[i][p] / A[p][p]
            for j in range(p + 1, n):
                A[i][j] -= f * A[p][j]
    return piv


def Wmat(s):
    s0, sU = s[0], s[1:]
    tot = s0 + sum(sU)
    W = [[Fr(0)] * (u + 1) for _ in range(u + 1)]
    W[0][0] = Fr(1)
    for i in range(u):
        for j in range(u):
            W[1 + i][1 + j] = s0 * ((sU[i] if i == j else 0) - sU[i] * sU[j] / tot)
    return W


ok = True
t0 = time.time()
pt = [Fr(1, 2), Fr(1, 4), Fr(1, 6), Fr(1, 12)]
piv = ldl(Gz(pt))
print(f'(P) G_z at s = {[str(x) for x in pt]}: LDL pivots {[float(x) for x in piv]}', flush=True)
ok &= len(piv) == u + 1 and all(x > 0 for x in piv)
# (C1)
h = Fr(1, 10)
found = None
for e in range(2, 12):
    for split in ((Fr(1, 3), Fr(1, 3), Fr(1, 3)), (Fr(8, 10), Fr(1, 10), Fr(1, 10)), (Fr(1, 2), Fr(1, 2 * 10 ** 3), Fr(1, 2) - Fr(1, 2 * 10 ** 3))):
        S = Fr(1, 10 ** e)
        s = [Fr(1)] + [S * x for x in split]
        G = Gz(s)
        W = Wmat(s)
        piv = ldl([[G[i][j] - h * W[i][j] for j in range(u + 1)] for i in range(u + 1)])
        if min(piv) <= 0:
            found = (s, piv)
            break
    if found:
        break
print(f'(C1) h = 1/10: G_z - h W not positive definite at s = {[str(x) for x in found[0]] if found else None}', flush=True)
ok &= found is not None
# (C2)
found2 = None
for e in range(1, 12):
    s = [Fr(1), Fr(1, 10 ** e), Fr(1, 10 ** e), Fr(1, 10 ** e)]
    G = Gz(s)
    G[0][0] -= Fr(1, 1000)
    piv = ldl(G)
    if min(piv) <= 0:
        found2 = s
        break
print(f'(C2) G_z - (1/1000) e_rho e_rho^T not positive definite at s = {[str(x) for x in found2] if found2 else None}', flush=True)
ok &= found2 is not None
print(f'time {time.time() - t0:.1f} s')
print('PASS' if ok else 'FAIL')
