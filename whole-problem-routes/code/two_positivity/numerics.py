"""Floating-point evidence for the remarks after Proposition 7.11 in whole-problem-routes/README.md (no claim rests on it)."""
import numpy as np
from common import *


def sect(Z):
    """Split an operator Z on (A1 A2) -> (B1 B2), viewed as a 9 x 9 matrix, into the four parts that are symmetric (S)
    or antisymmetric (A) under the swap A_k <-> B_k at copy 1 and at copy 2."""
    T = Z.reshape(3, 3, 3, 3)  # a1 a2 b1 b2
    ZT1 = T.transpose(2, 1, 0, 3).reshape(9, 9)
    ZT2 = T.transpose(0, 3, 2, 1).reshape(9, 9)
    ZT12 = T.transpose(2, 3, 0, 1).reshape(9, 9)
    SS = (Z + ZT1 + ZT2 + ZT12) / 4; SA = (Z + ZT1 - ZT2 - ZT12) / 4
    AS = (Z - ZT1 + ZT2 - ZT12) / 4; AA = (Z - ZT1 - ZT2 + ZT12) / 4
    return SS, SA, AS, AA

rng = np.random.default_rng(2026)
C = choi(0.75)
W = two_copy(C, C)
# (N1) SR2 minimum of C(x)C at a = 3/4, many complex and real starts
m1 = sr2_min(W, D=9, starts=150, rng=rng)
m2 = sr2_min(W, D=9, starts=50, rng=rng, real=True)
print(f'(N1) a=3/4: min over SR2 of <C(x)C>: complex starts {m1:.3e}, real starts {m2:.3e}')
# (N2) X_w (x) X_w: SR2-positive exactly for w in [1, 2]
def X(w): return M + (w - 1) * P_an
for w in (0.95, 0.99, 1.0, 1.5, 2.0, 2.01, 2.05):
    print(f'(N2) w={w}: min SR2 of X_w(x)X_w = {sr2_min(two_copy(X(w), X(w)), D=9, starts=24, rng=rng):.3e}')
# (N3) zero family at the w -> 2+ limit: minimiser at w = 2.001, evaluated for X_2 (x) X_2 = 4 C(x)C
m, Z = sr2_min(two_copy(X(2.001), X(2.001)), D=9, starts=24, rng=rng, return_vec=True)
Zn = Z / np.linalg.norm(Z)
print(f'(N3) w=2.001 min {m:.3e}; value of 4C(x)C there {form(4 * W, Zn):.3e}; sectors SS,SA,AS,AA =',
      [round(np.linalg.norm(S) ** 2, 4) for S in sect(Zn)])
U, s, Vh = np.linalg.svd(Zn); u = U[:, :2] * s[:2]; v = Vh[:2].T
for vec, nm in ((u, 'A'), (v, 'B')):
    t = vec.reshape(3, 3, 2)
    e1 = np.linalg.eigvalsh(np.einsum('ijr,kjr->ik', t, t.conj())); e2 = np.linalg.eigvalsh(np.einsum('ijr,ikr->jk', t, t.conj()))
    print(f'     {nm}1 marginal eigenvalues {np.round(e1, 4)}, {nm}2 marginal eigenvalues {np.round(e2, 4)}')
# (N4) product zero Q (x) x x^T
Q = np.diag([1., 1., 0.]); x = np.array([0.6, 0.0, 0.8]); Z0 = np.kron(Q, np.outer(x, x))
print(f'(N4) <C(x)C> at Q(x)xx^T (normalised): {form(W, Z0):.3e}; rank of Z0: {np.linalg.matrix_rank(Z0)}')
# (N5) restrictions of X_w to a plane (random complex isometry V: C^2 -> C^3, on the A side) are PSD and PPT, w in [1,2]
worst = np.inf
for trial in range(200):
    V = np.linalg.qr(rng.normal(size=(3, 2)) + 1j * rng.normal(size=(3, 2)))[0]
    L = np.kron(V, np.eye(3))
    for w in (1.0, 1.5, 2.0):
        R = L.conj().T @ X(w) @ L                       # on C^2 (x) C^3
        RG = R.reshape(2, 3, 2, 3).transpose(0, 3, 2, 1).reshape(6, 6)   # partial transpose on the C^3 factor
        worst = min(worst, np.linalg.eigvalsh(R).min(), np.linalg.eigvalsh(RG).min())
print(f'(N5) min eigenvalue of X_w|plane and of its partial transpose over 200 random complex planes, w in {{1,1.5,2}}: {worst:.3e}')
