"""Shared helpers for the binary-merger (BM) tests.

rho = I - F/2 on C^3 (x) C^3,  G(A, B) = (A (x) B) rho (A (x) B)^*  (a "filtered Werner atom").
K_1(3) = cone{G(A, B)} (products are atoms).  SF: every element of K_1(3) is lambda G(E, F) + separable.
BM: G_1 + G_2 = G_3 + separable for every two atoms.

Operators on C^3 (x) C^3 are 9 x 9 arrays, index 3*i + j (first factor i, second factor j).
"""
import numpy as np
import torch
from scipy.optimize import nnls

torch.set_default_dtype(torch.float64)
D = 3
N = 9
F9 = np.zeros((N, N))
for i in range(D):
    for j in range(D):
        F9[3 * i + j, 3 * j + i] = 1.0
RHO = np.eye(N) - F9 / 2
RHO_T = torch.tensor(RHO, dtype=torch.complex128)
OMEGA = np.eye(D).reshape(N)            # sum_i |ii>


def G(A, B):
    K = np.kron(A, B)
    Y = K @ RHO @ K.conj().T
    return (Y + Y.conj().T) / 2


def Gt(A, B):
    K = torch.kron(A, B)
    return K @ RHO_T @ K.conj().T


def pt(Y):
    """partial transpose on the second factor"""
    if isinstance(Y, np.ndarray):
        return Y.reshape(D, D, D, D).transpose(0, 3, 2, 1).reshape(N, N)
    return Y.reshape(D, D, D, D).permute(0, 3, 2, 1).reshape(N, N)


def ptmin(Y):
    Y = (Y + Y.conj().T) / 2
    return np.linalg.eigvalsh(pt(Y))[0] / np.trace(Y).real


IU = np.triu_indices(N, 1)


def vecr(M):
    """Frobenius-isometric real vectorisation of a Hermitian 9 x 9 matrix (sqrt 2 on off-diagonal entries)"""
    return np.concatenate([np.diag(M).real, np.sqrt(2) * M[IU].real, np.sqrt(2) * M[IU].imag])


def vecr_batch(Ms):
    """Ms: (K, 9, 9) -> (81, K)"""
    d = np.real(np.einsum('kii->ki', Ms))
    up = Ms[:, IU[0], IU[1]]
    return np.concatenate([d, np.sqrt(2) * up.real, np.sqrt(2) * up.imag], axis=1).T


def unvecr(v):
    M = np.zeros((N, N), complex)
    M[np.diag_indices(N)] = v[:N]
    m = len(IU[0])
    M[IU] = (v[N:N + m] + 1j * v[N + m:]) / np.sqrt(2)
    M = M + np.triu(M, 1).conj().T
    return M


def rand_c(g, *shape):
    return g.standard_normal(shape) + 1j * g.standard_normal(shape)


def product_cols(xs, ys):
    """xs, ys: (K, 3) -> normalised product projectors (K, 9, 9)"""
    xs = xs / np.linalg.norm(xs, axis=1, keepdims=True)
    ys = ys / np.linalg.norm(ys, axis=1, keepdims=True)
    v = np.einsum('ka,kb->kab', xs, ys).reshape(-1, N)
    return np.einsum('ka,kb->kab', v, v.conj())


def best_products(R, g, nstart=24, iters=60, keep=8):
    """local maximisers of <xy|R|xy> over unit x, y (alternating top eigenvectors); returns list (value, x, y)"""
    R4 = R.reshape(D, D, D, D)
    xs = rand_c(g, nstart, D); ys = rand_c(g, nstart, D)
    ys /= np.linalg.norm(ys, axis=1, keepdims=True)
    for _ in range(iters):
        Mx = np.einsum('kb,abcd,kd->kac', ys.conj(), R4, ys)
        w, U = np.linalg.eigh(Mx); xs = U[:, :, -1]
        My = np.einsum('ka,abcd,kc->kbd', xs.conj(), R4, xs)
        w, U = np.linalg.eigh(My); ys = U[:, :, -1]
    vals = w[:, -1]
    order = np.argsort(-vals)
    out = []
    for k in order[:keep]:
        out.append((vals[k], xs[k].copy(), ys[k].copy()))
    return out


def min_product_value(R, g, nstart=24, iters=80):
    """min over unit product vectors of <xy|R|xy> (local searches); for separability-witness checks"""
    vals = best_products(-R, g, nstart, iters, keep=1)
    return -vals[0][0]


class SepDict:
    """dictionary of normalised product projectors in vecr form (the product vectors are kept in xv, yv)"""

    def __init__(self, g, nrand=600):
        self.g = g
        xs, ys = rand_c(g, nrand, D), rand_c(g, nrand, D)
        self.xv = xs / np.linalg.norm(xs, axis=1, keepdims=True)
        self.yv = ys / np.linalg.norm(ys, axis=1, keepdims=True)
        self.V = vecr_batch(product_cols(self.xv, self.yv))

    def add(self, x, y):
        self.add_many([x], [y])

    def add_many(self, xs, ys):
        xs = np.asarray(xs); ys = np.asarray(ys)
        xs = xs / np.linalg.norm(xs, axis=1, keepdims=True); ys = ys / np.linalg.norm(ys, axis=1, keepdims=True)
        self.xv = np.concatenate([self.xv, xs]); self.yv = np.concatenate([self.yv, ys])
        self.V = np.concatenate([self.V, vecr_batch(product_cols(xs, ys))], axis=1)

    def prune(self, coef, keep_extra=300):
        """keep active columns plus a random subset"""
        act = np.where(coef > 0)[0]
        rest = np.setdiff1d(np.arange(self.V.shape[1]), act)
        if len(rest) > keep_extra:
            rest = self.g.choice(rest, keep_extra, replace=False)
        keep = np.concatenate([act, rest])
        self.V = self.V[:, keep]; self.xv = self.xv[keep]; self.yv = self.yv[keep]


def nnls_fit(target_vec, cols):
    coef, rn = nnls(cols, target_vec, maxiter=50 * cols.shape[1])
    return coef, rn
