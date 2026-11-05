"""Checks for Section 6 of NOTES.md: the three-copy Kronecker example (stacked ranks of A = span{g, h} and
B_eps = span{w1 + eps g, w2 + eps h} at every copy, the 6 x 6 minor = eps^3 exactly, g an anchor, B_eps anchor-free,
A balanced), and the range compressions of the GHZ plane and of a generic local-qubit plane (unique zero; next
eigenvalue 2^{1-n} for GHZ), n = 2, 3."""
import numpy as np, itertools, sys
from fractions import Fraction as Fr
import sympy as sp
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_anchors_matrixfree import build, min_eig_complement, M1, apply_local
def idx(s):
    v = 0
    for c in s: v = 3 * v + c
    return v
n = 3; N = 27
def vec(terms):
    v = np.zeros(N)
    for s, c in terms: v[idx(s)] += c
    return v
w1 = vec([((1, 0, 0), 1), ((0, 1, 0), 1), ((0, 0, 1), 1)]); w2 = vec([((2, 0, 0), 1), ((0, 2, 0), 1), ((0, 0, 2), 1)])
g = vec([((0, 0, 0), 1), ((1, 1, 1), 1), ((2, 2, 2), 1)]); h = vec([((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1)])
eps = 1e-8
# (i) stacked rank at each copy, exact with eps = 1/10^8
def stacked(vs, copy, e):
    rows = []
    for v in vs:
        T = np.array(v, dtype=object).reshape(3, 3, 3); T = np.moveaxis(T, copy, 0).reshape(3, 9)
        rows += [list(T[r]) for r in range(3)]
    return sp.Matrix(rows)
E = sp.Rational(1, 10 ** 8)
gs = [sp.Integer(int(x)) for x in g]; hs = [sp.Integer(int(x)) for x in h]
w1s = [sp.Integer(int(x)) for x in w1]; w2s = [sp.Integer(int(x)) for x in w2]
Bv = [[w1s[i] + E * gs[i] for i in range(N)], [w2s[i] + E * hs[i] for i in range(N)]]
Av = [gs, hs]
ranks = [(stacked(Av, c, E).rank(), stacked(Bv, c, E).rank()) for c in range(3)]
# the 6 x 6 minor at column words 00,01,02,10,20,22 for the first copy
S = stacked(Bv, 0, E); colsel = [idx2 for idx2 in [0, 1, 2, 3, 6, 8]]
minor = S[:, colsel].det()
print(f"V5: stacked ranks (A, B_eps) at the three copies: {ranks}; the 6x6 minor det = {minor} (eps^3 = {E**3})")
# (ii) anchors: g strict anchor (Theorem 12), B_eps anchor-free
rng = np.random.default_rng(0); worst = -9
for trial in range(4):
    c = rng.normal(size=2) + 1j * rng.normal(size=2); p = c[0] * (w1 + eps * g) + c[1] * (w2 + eps * h); p = p / np.linalg.norm(p)
    lm = min_eig_complement(build(p, n), tol=1e-10); worst = max(worst, lm)
lg = min_eig_complement(build(g / np.linalg.norm(g), n), tol=1e-10)
print(f"    max over 4 unit vectors of B_eps of lambda_min(Y_p) off kernel = {worst:+.5f} (claim < -2/387 = {-2/387:+.5f}); g: {lg:+.5f} (anchor)")
# (iii) marginals of A are I/3 for random unit vectors
c = rng.normal(size=2) + 1j * rng.normal(size=2); p = (c[0] * g + c[1] * h); p /= np.linalg.norm(p)
T = p.reshape(3, 3, 3); marg = max(np.abs(np.einsum('aij,bij->ab', np.moveaxis(T, k, 0), np.moveaxis(T, k, 0).conj()) - np.eye(3) / 3).max() for k in range(3))
print(f"    A: one-copy marginals of a random unit vector, max |rho - I/3| = {marg:.1e}")
# ---------- GHZ plane compression and generic local-qubit plane compression (n = 2, 3) ----------
def compression(U, nn):
    NN = 3 ** nn; m = M1(); Q = np.zeros((2 * NN, 2 * NN), complex)
    for j in range(2):
        for xc in range(NN):
            Z = np.zeros((NN, NN), complex); Z[:, xc] = U[:, j]
            MZ = apply_local(Z, m, nn)
            Q[:, j * NN + xc] = (U.conj().T @ MZ).reshape(-1)
    return np.linalg.eigvalsh((Q + Q.conj().T) / 2)
for nn in (2, 3):
    NN = 3 ** nn
    U = np.zeros((NN, 2)); U[0, 0] = 1; U[NN - 1 - (NN - 1) // 2 * 0, 1] = 0
    U = np.zeros((NN, 2)); U[0, 0] = 1; U[sum(3 ** i for i in range(nn)), 1] = 1     # |0..0>, |1..1>
    ev = compression(U, nn)
    print(f"GHZ plane n={nn}: compression eigenvalues: lowest {ev[:3]}, #zero = {np.sum(np.abs(ev) < 1e-12)}; gamma_n = {2.0**(1-nn)}")
    # generic local-qubit plane
    qub = [i for i, s in enumerate(itertools.product(range(3), repeat=nn)) if 2 not in s]
    A = np.zeros((NN, 2), complex); A[qub] = rng.normal(size=(len(qub), 2)) + 1j * rng.normal(size=(len(qub), 2))
    U2, _ = np.linalg.qr(A); ev2 = compression(U2, nn)
    print(f"generic local-qubit plane n={nn}: lowest compression eigenvalues {ev2[:3]}, #zero = {np.sum(np.abs(ev2) < 1e-10)}")
