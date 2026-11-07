"""Shared floating-point utilities (exploration only; no claim rests on these).

Conventions.  One copy: C^3_A (x) C^3_B, index 3*a + b.  A vector chi in A (x) B is vec(Z), Z_{ab}.
Two copies: chi = vec(Z), Z a 9 x 9 matrix with rows (a1, a2) and columns (b1, b2); vector index
((a1*3 + a2)*3 + b1)*3 + b2.
"""
import numpy as np

d = 3
I3 = np.eye(3)
I9 = np.eye(9)
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3):
        F[3 * i + j, 3 * j + i] = 1.0
Om = np.zeros(9)
for i in range(3):
    Om[4 * i] = 1.0
OO = np.outer(Om, Om)
P_Om = OO / 3
P_sym = (I9 + F) / 2
P_an = (I9 - F) / 2
P_s0 = P_sym - P_Om


def choi(a):
    return a * I9 - (F + OO) / 4


M = I9 - OO / 2


def two_copy(W1, W2):
    """W1 on A1B1, W2 on A2B2  ->  operator on A1 A2 B1 B2 (index order a1, a2, b1, b2)."""
    T = np.kron(W1, W2).reshape([3] * 8)  # (a1 b1 a2 b2 ; a1' b1' a2' b2')
    T = T.transpose(0, 2, 1, 3, 4, 6, 5, 7)
    return T.reshape(81, 81)


def sr2_min(W, D=9, starts=10, iters=500, rng=None, real=False, tol=1e-14, return_vec=False):
    """min of <chi|W|chi> over unit chi of Schmidt rank <= 2 in C^D (x) C^D (W is D^2 x D^2), alternating."""
    if rng is None:
        rng = np.random.default_rng(0)
    W4 = W.reshape(D, D, D, D)  # (a b ; a' b')
    best = np.inf
    bestZ = None
    for st in range(starts):
        if real:
            U = rng.normal(size=(D, 2))
            V = rng.normal(size=(D, 2))
        else:
            U = rng.normal(size=(D, 2)) + 1j * rng.normal(size=(D, 2))
            V = rng.normal(size=(D, 2)) + 1j * rng.normal(size=(D, 2))
        val = np.inf
        for it in range(iters):
            # fix U (orthonormal), optimise V: chi = sum_r u_r (x) v_r
            U, _ = np.linalg.qr(U)
            # H[(r,b),(s,b')] = sum_{a a'} conj(U[a,r]) W[a b, a' b'] U[a', s]
            H = np.einsum('ar,abcd,cs->rbsd', U.conj(), W4, U).reshape(2 * D, 2 * D)
            H = (H + H.conj().T) / 2
            w, X = np.linalg.eigh(H)
            V = X[:, 0].reshape(2, D).T
            V, _ = np.linalg.qr(V)
            H = np.einsum('br,abcd,ds->rasc', V.conj(), W4, V).reshape(2 * D, 2 * D)
            H = (H + H.conj().T) / 2
            w2, X2 = np.linalg.eigh(H)
            U = X2[:, 0].reshape(2, D).T
            if abs(w2[0] - val) < tol:
                val = w2[0]
                break
            val = w2[0]
        if val < best:
            best = val
            bestZ = U @ V.T
    if return_vec:
        return best, bestZ
    return best


def form(W, Z):
    z = Z.reshape(-1)
    return np.real(z.conj() @ W @ z) / np.real(z.conj() @ z)
