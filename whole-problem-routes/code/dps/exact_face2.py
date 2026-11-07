"""Exact certificate at t = 4^-n (facial reduction), Gamma_A-only form -- faster version of exact_face.py.

Same mathematics as exact_face.py: unknowns are orbit coefficients u of Q_A_raw (orbits of diagram tuples under copy
permutations x transpose); equations G_beta(Q_A) f_beta = 0 and G_beta(Z) f_beta = 0 for the forced vectors f_beta,
beta in {1,3}^n.  Rows are computed through adjoints of the coordinate permutations acting on product tensors (int64,
with an explicit overflow bound), independent rows are selected modulo a prime, and the minimal correction is solved
exactly over Q (python-flint); all equations are then checked exactly.

usage: python exact_face2.py sol.npz n t_num t_den log2L out.txt     (a = 3/4)
"""
import sys
import math
import itertools
from fractions import Fraction
import numpy as np
import sympy as sp
import flint
import brauer as br

sol, n = sys.argv[1], int(sys.argv[2])
t = Fraction(int(sys.argv[3]), int(sys.argv[4]))
log2L = int(sys.argv[5])
out = sys.argv[6]
NB = br.NDIAG
d = np.load(sol)
nm = br.named()
piA = np.array(br.pt_perm([0]))
tr = np.array(br.pt_perm([0, 1, 2]))
muL = np.array([br.MULT[nm['S_BB']][x][0] for x in range(NB)])
muR = np.array([br.MULT[x][nm['S_BB']][0] for x in range(NB)])
muRL = muR[muL]

allk = np.array(list(itertools.product(range(NB), repeat=n)), dtype=np.int64)
pw = NB ** np.arange(n - 1, -1, -1)
canon = np.minimum(np.sort(allk, axis=1) @ pw, np.sort(tr[allk], axis=1) @ pw)
orb_ids, orb = np.unique(canon, return_inverse=True)
NO = len(orb_ids)
print('orbits', NO, flush=True)

Bex = br.spin_m0_bases()
B = {j: np.array(Bex[j].tolist(), dtype=np.int64) for j in range(4)}
Psi = {j: np.array([B[j].T @ br.DMAT[a] @ B[j] for a in range(NB)]) for j in range(4)}
Ssym = np.zeros((27, 27), dtype=np.int64)
for p in itertools.permutations(range(3)):
    for i in itertools.product(range(3), repeat=3):
        j_ = tuple(i[p[k]] for k in range(3))
        Ssym[9 * j_[0] + 3 * j_[1] + j_[2], 9 * i[0] + 3 * i[1] + i[2]] += 1
fhat = {}
for j in (1, 3):
    Bm = sp.Matrix(B[j].tolist())
    img = sp.Matrix((Ssym @ B[j]).tolist())
    assert img.rank() == 1
    col = [k for k in range(img.shape[1]) if any(img[:, k])][0]
    s_ = Bm.solve_least_squares(img[:, col])
    assert Bm * s_ == img[:, col]
    s_ = s_ * sp.ilcm(1, 1, *[sp.fraction(v)[1] for v in s_])
    g = sp.igcd(0, 0, *[int(v) for v in s_]) if len(s_) > 1 else abs(int(s_[0]))
    fhat[j] = np.array([int(v) // g for v in s_], dtype=np.int64)
V = {j: np.array([Psi[j][a] @ fhat[j] for a in range(NB)], dtype=np.int64) for j in (1, 3)}   # (15, m_j)
vmax = max(int(np.abs(V[j]).max()) for j in (1, 3))
omax = int(np.bincount(orb).max())
bound = vmax ** n * 16 * omax
print('fhat', {j: fhat[j].tolist() for j in fhat}, 'max|V|', vmax, 'overflow bound', bound, '< 2^63:', bound < 2 ** 63,
      flush=True)
assert bound < 2 ** 63


def adj_perm(K, pi):
    """adjoint of Y[pi(a)] = X[a]:  (K o pi)[a] = K[pi(a1), .., pi(an)]"""
    return K[np.ix_(*([pi] * n))]


def adj_Pc4(K):
    return K + adj_perm(K, muL) + adj_perm(K, muR) + adj_perm(K, muRL)


def orbit_sum(T):
    acc = np.zeros(NO, dtype=np.int64)
    np.add.at(acc, orb, T.reshape(-1))
    return acc


cvec = [Fraction(0)] * NB
cvec[nm['I']] += Fraction(3, 4)
cvec[nm['F_AB']] -= Fraction(1, 4)
cvec[nm['Om_AB']] -= Fraction(1, 4)


def x_dot_K1(beta, r):
    """<X, K1> exactly, X = W(x)I - t I, K1 = adj_Pc4(K): a sum of four product tensors"""
    tot = Fraction(0)
    for g in (np.arange(NB), muL, muR, muRL):
        pw_, pi_ = Fraction(1), 1
        for s in range(n):
            w = V[beta[s]][g, r[s]]                       # w[b] = V[g(b), r_s]
            pw_ *= sum(cvec[b] * int(w[b]) for b in range(NB))
            pi_ *= int(w[nm['I']])
        tot += pw_ - t * pi_
    return tot


rows1, rows2, rhs2 = [], [], []
for beta in itertools.product((1, 3), repeat=n):
    for r in itertools.product(*[range(V[j].shape[1]) for j in beta]):
        K = np.array([1], dtype=np.int64)
        for s in range(n):
            K = np.multiply.outer(K, V[beta[s]][:, r[s]]).reshape(-1)
        K = K.reshape([NB] * n)
        K1 = adj_Pc4(K)                         # functional of 4 P(.)P
        rows1.append(orbit_sum(K1))             # = 4 G(Q_A) f   (on x = E u)
        K2 = adj_Pc4(adj_perm(K1, piA))         # functional of 16 P Gamma_A(P(.)P) P
        rows2.append(orbit_sum(K2))
        rhs2.append(4 * x_dot_K1(beta, r))      # 16 G(P X P) f = 4 * <X, K1>
A = np.array(rows1 + rows2, dtype=object)
bvec = [Fraction(0)] * len(rows1) + [Fraction(v) for v in rhs2]
print('equations', A.shape, flush=True)
# independent rows modulo a prime
p = 2 ** 61 - 1
Amod = flint.nmod_mat(A.shape[0], A.shape[1], [int(v) % p for v in A.reshape(-1)], p)
R, rk = Amod.transpose().rref()
piv = []
for i in range(rk):
    piv.append(next(j for j in range(A.shape[0]) if int(R[i, j]) != 0))
print('rank (mod p)', rk, flush=True)
Ar = flint.fmpq_mat(len(piv), NO, [int(v) for v in A[piv].reshape(-1)])
br_ = flint.fmpq_mat(len(piv), 1, [flint.fmpq(bvec[i].numerator, bvec[i].denominator) for i in piv])
cnt = np.bincount(orb)
uf = np.bincount(orb, weights=d['xA']) / cnt
L = 2 ** log2L
u0 = flint.fmpq_mat(NO, 1, [flint.fmpq(int(round(v * L)), L) for v in uf])
res = br_ - Ar * u0
G = Ar * Ar.transpose()
mu = G.solve(res)
u = u0 + Ar.transpose() * mu
# check ALL equations exactly (in chunks)
for c0 in range(0, A.shape[0], 256):
    rows_ = list(range(c0, min(A.shape[0], c0 + 256)))
    Ach = flint.fmpq_mat(len(rows_), NO, [int(v) for v in A[rows_].reshape(-1)])
    bch = flint.fmpq_mat(len(rows_), 1, [flint.fmpq(bvec[i].numerator, bvec[i].denominator) for i in rows_])
    assert Ach * u == bch, 'not all forced-kernel equations hold'
print('all', A.shape[0], 'equations hold exactly; max correction',
      max(abs(float(u[i, 0] - u0[i, 0])) for i in range(NO)), flush=True)
uq = [Fraction(int(u[i, 0].p), int(u[i, 0].q)) for i in range(NO)]
D = 1
for v in uq:
    D = math.lcm(D, v.denominator)
XA = np.array([int(v * D) for v in uq], dtype=object)[orb]
with open(out, 'w') as f:
    f.write(f'n {n}\n')
    f.write('a 3 4\n')
    f.write(f't {t.numerator} {t.denominator}\n')
    f.write(f'L {D}\n')
    for k, m in enumerate(br.DIAGRAMS):
        f.write(f'D{k} ' + ' '.join(f'{p_} {q}' for p_, q in m) + '\n')
    f.write('QA ' + ' '.join(map(str, XA.tolist())) + '\n')
    f.write('QB ' + ' '.join(['0'] * len(XA)) + '\n')
print('wrote', out, 'denominator bits', D.bit_length())
