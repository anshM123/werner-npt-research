"""Exact certificate AT the product value t = 4^-n (facial reduction), Gamma_A-only form (Q_B = 0).

At t = 4^-n every real product minimiser x = y = v_1 (x) .. (x) v_n forces Q_A and Z to annihilate
psi = x (x) y (x) y; these span, per copy, the totally symmetric tensors Sym^3(C^3), whose m = 0 vectors are the whole
spin-3 space and one vector of the spin-1 space.  So in each block beta in {1,3}^n, Q_A and Z must kill one forced
vector f_beta.  This script takes a floating solution (sym_sdp_cones.py ... A), averages it over the symmetry orbits
(copy permutations x transpose), rounds, and then moves it by an exact minimal correction onto the rational affine
space {G_beta(Q_A) f_beta = 0, G_beta(Z) f_beta = 0 for all beta in {1,3}^n}, where Z = P[W(x)I - t I - Gamma_A(Q_A)]P.
The result is written as a certificate (QB = 0) for exact_verify.py, whose LDL^T accepts the exact zero pivots.

usage: python exact_face.py sol.npz n t_num t_den log2L out.txt     (a = 3/4)
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
for x in range(NB):
    assert br.MULT[nm['S_BB']][x][1] == 0 and br.MULT[x][nm['S_BB']][1] == 0


def perm(X, pi):
    Y = np.empty_like(X)
    Y[np.ix_(*([pi] * n))] = X
    return Y


def Pc4(X):            # 4 P X P
    SX = perm(X, muL)
    return X + SX + perm(X, muR) + perm(SX, muR)


# orbits of diagram tuples under copy permutations x transpose
allk = np.array(list(itertools.product(range(NB), repeat=n)), dtype=np.int64)
pw = NB ** np.arange(n - 1, -1, -1)
canon = np.minimum(np.sort(allk, axis=1) @ pw, np.sort(tr[allk], axis=1) @ pw)
orb_ids, orb = np.unique(canon, return_inverse=True)
NO = len(orb_ids)
print('orbits', NO)

# per-copy exact bases and forced vectors
Bex = br.spin_m0_bases()
B = {j: np.array(Bex[j].tolist(), dtype=object) for j in range(4)}
Psi = {j: [B[j].T.dot(br.DMAT[a].astype(object)).dot(B[j]) for a in range(NB)] for j in range(4)}
# symmetriser over the 3 strands of one copy
Ssym = np.zeros((27, 27), dtype=object)
for p in itertools.permutations(range(3)):
    for i in itertools.product(range(3), repeat=3):
        j_ = tuple(i[p[k]] for k in range(3))
        Ssym[9 * j_[0] + 3 * j_[1] + j_[2], 9 * i[0] + 3 * i[1] + i[2]] += 1
fhat = {}
for j in (1, 3):
    img = sp.Matrix(Ssym.dot(B[j]).tolist())
    assert img.rank() == 1
    col = [k for k in range(img.shape[1]) if any(img[:, k])][0]
    f = img[:, col]
    sol_ = sp.Matrix(B[j].tolist()).solve_least_squares(f)
    assert sp.Matrix(B[j].tolist()) * sol_ == f
    sol_ = sol_ * sp.ilcm(1, 1, *[sp.fraction(v)[1] for v in sol_])
    fhat[j] = [int(v) for v in sol_]
    print('forced m=0 vector, spin', j, 'coordinates', fhat[j])
# v_j(a) = Psi_j(a) fhat_j
V = {j: np.array([Psi[j][a].dot(np.array(fhat[j], dtype=object)) for a in range(NB)], dtype=object) for j in (1, 3)}
forced_blocks = list(itertools.product((1, 3), repeat=n))


def functionals(Y):
    """[G_beta(Y) fhat_beta]_r for all forced beta and r, as a list of integers (Y integer tensor)"""
    outv = []
    for beta in forced_blocks:
        cur = Y
        for s in reversed(range(n)):
            cur = np.tensordot(cur, V[beta[s]], axes=([s], [0]))     # contracts a_s, appends r_s
        outv += list(np.transpose(cur, list(range(n - 1, -1, -1))).reshape(-1))
    return outv


# constraint rows on orbit coordinates u (x = E u): Q_A = P x P,  Z = P (W(x)I - t - Gamma_A(Q_A)) P
rows1, rows2 = [], []
for o in range(NO):
    e = (orb == o).astype(np.int64).astype(object).reshape([NB] * n)
    qa4 = Pc4(e)                               # 4 Q_A
    rows1.append(functionals(qa4))             # G(4 Q_A) f
    z16 = Pc4(perm(qa4, piA))                  # 16 P Gamma_A(Q_A) P
    rows2.append(functionals(z16))
A1 = np.array(rows1, dtype=object).T           # (#func, NO):  A1 u = 4 G(Q_A) f
A2 = np.array(rows2, dtype=object).T           # A2 u = 16 G(P Gamma_A(Q_A) P) f
# constant part: G(P(W(x)I - tI)P) f
cvec = [Fraction(0)] * NB
cvec[nm['I']] += Fraction(3, 4)
cvec[nm['F_AB']] -= Fraction(1, 4)
cvec[nm['Om_AB']] -= Fraction(1, 4)
Xc = np.array([Fraction(1)], dtype=object)
for _ in range(n):
    Xc = np.multiply.outer(Xc, np.array(cvec, dtype=object)).reshape(-1)
Xc = Xc.reshape([NB] * n)
Ic = np.zeros([NB] * n, dtype=object)
Ic[tuple([nm['I']] * n)] = 1
Xc = Xc - t * Ic
den = 1
for v in Xc.reshape(-1):
    den = math.lcm(den, Fraction(v).denominator)
Xi = np.array([int(Fraction(v) * den) for v in Xc.reshape(-1)], dtype=object).reshape([NB] * n)
bX = functionals(Pc4(Xi))                      # = 4 den G(X) f
# equations:  A1 u = 0 ;  (16 G(X) f) - A2 u = 0  ->  A2 u = 16 G(X) f = 4 * bX / den
Aeq = np.vstack([A1, A2])
beq = [Fraction(0)] * A1.shape[0] + [Fraction(4 * v, den) for v in bX]
print('equations', Aeq.shape)
# rounded orbit-averaged float solution
cnt = np.bincount(orb)
uf = np.bincount(orb, weights=d['xA']) / cnt
L = 2 ** log2L
u0 = [Fraction(int(round(v * L)), L) for v in uf]
Am = flint.fmpq_mat([[int(v) for v in row] for row in Aeq])
u0m = flint.fmpq_mat([[v.numerator, ] for v in u0]) if False else flint.fmpq_mat(NO, 1, [flint.fmpq(v.numerator, v.denominator) for v in u0])
bm = flint.fmpq_mat(len(beq), 1, [flint.fmpq(v.numerator, v.denominator) for v in beq])
r = bm - Am * u0m
print('residual before correction (max abs):', max(abs(float(r[i, 0])) for i in range(r.nrows())))
G = Am * Am.transpose()
aug = flint.fmpq_mat(G.nrows(), G.ncols() + 1, [G[i, j] if j < G.ncols() else r[i, 0]
                                                for i in range(G.nrows()) for j in range(G.ncols() + 1)])
R, rk = aug.rref()
print('rank of normal system', rk)
mu = [flint.fmpq(0)] * G.ncols()
for i in range(rk):
    piv = next(j for j in range(G.ncols() + 1) if R[i, j] != 0)
    assert piv < G.ncols(), 'inconsistent system: the forced-kernel equations have no solution'
    mu[piv] = R[i, G.ncols()]
mum = flint.fmpq_mat(G.ncols(), 1, mu)
u = u0m + Am.transpose() * mum
assert Am * u == bm
corr = max(abs(float(u[i, 0] - u0m[i, 0])) for i in range(NO))
print('exact solution found; max correction', corr)
# common denominator and output
uf_ = [Fraction(int(u[i, 0].p), int(u[i, 0].q)) for i in range(NO)]
D = 1
for v in uf_:
    D = math.lcm(D, v.denominator)
U_int = [int(v * D) for v in uf_]
XA = np.array(U_int, dtype=object)[orb]
with open(out, 'w') as f:
    f.write(f'n {n}\n')
    f.write('a 3 4\n')
    f.write(f't {t.numerator} {t.denominator}\n')
    f.write(f'L {D}\n')
    for k, m in enumerate(br.DIAGRAMS):
        f.write(f'D{k} ' + ' '.join(f'{p} {q}' for p, q in m) + '\n')
    f.write('QA ' + ' '.join(map(str, XA.tolist())) + '\n')
    f.write('QB ' + ' '.join(['0'] * len(XA)) + '\n')
print('wrote', out, 'denominator bits', D.bit_length())
