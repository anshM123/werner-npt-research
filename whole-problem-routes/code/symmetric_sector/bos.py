"""Bosonic (permutation-symmetric) sector of P_a^{(x)n}, P_a(X) = a Tr(X) I - (X + X^T)/4 on M_3.

K_n = (Pi (x) Pi) C_a^{(x)n} (Pi (x) Pi) on Sym^n(R^3) (x) Sym^n(R^3), in the orthonormal monomial basis
e_alpha (|alpha| = n), built by the recursion through the isometric embedding
    Sym^n -> C^3 (x) Sym^{n-1},  e_alpha -> sum_i e_i (x) L_i e_alpha,   L_i e_alpha = sqrt(alpha_i/n) e_{alpha - d_i}.
Then K_n = sum_{ijkl} C_{ij,kl} (L_i (x) L_j)^T K_{n-1} (L_k (x) L_l),  C_{ij,kl} = a d_ik d_jl - d_il d_jk/4 - d_ij d_kl/4.
V(x, y) = (x (x) y)^T K_n (x (x) y) = <y| P_a^{(x)n}(x x^T) |y> for real symmetric x, y.
"""
import numpy as np
from functools import lru_cache
from math import comb, factorial, sqrt


@lru_cache(None)
def mindices(n, d=3):
    if d == 1:
        return ((n,),)
    out = []
    for k in range(n, -1, -1):
        for rest in mindices(n - k, d - 1):
            out.append((k,) + rest)
    return tuple(out)


@lru_cache(None)
def index_of(n, d=3):
    return {al: i for i, al in enumerate(mindices(n, d))}


def lowering(n, d=3):
    """For each i: (target index array into level n-1 basis (or -1), scale array sqrt(alpha_i/n))."""
    idx = index_of(n - 1, d)
    res = []
    for i in range(d):
        tgt = []
        sc = []
        for al in mindices(n, d):
            if al[i] > 0:
                be = list(al); be[i] -= 1
                tgt.append(idx[tuple(be)]); sc.append(sqrt(al[i] / n))
            else:
                tgt.append(0); sc.append(0.0)
        res.append((np.array(tgt), np.array(sc)))
    return res


def build_K(n, a, d=3):
    """K_n as a dense (N^2 x N^2) array, pair index (alpha, beta) -> alpha*N + beta."""
    K = np.array([[1.0]])  # n = 0
    for m in range(1, n + 1):
        Nm1 = len(mindices(m - 1, d)); N = len(mindices(m, d))
        low = lowering(m, d)
        # pair maps: r[i][j] (row index into K_{m-1}), s[i][j] (scale)
        r = {}; s = {}
        for i in range(d):
            ti, si = low[i]
            for j in range(d):
                tj, sj = low[j]
                r[i, j] = (ti[:, None] * Nm1 + tj[None, :]).reshape(-1)
                s[i, j] = (si[:, None] * sj[None, :]).reshape(-1)
        Kn = np.zeros((N * N, N * N))
        for i in range(d):
            for j in range(d):
                # identity term: a * (L_i L_j)^T K (L_i L_j)
                Kn += a * s[i, j][:, None] * s[i, j][None, :] * K[np.ix_(r[i, j], r[i, j])]
                # swap term: -(1/4) (L_i L_j)^T K (L_j L_i)
                Kn -= 0.25 * s[i, j][:, None] * s[j, i][None, :] * K[np.ix_(r[i, j], r[j, i])]
        # Omega term: -(1/4) (sum_i L_i L_i)^T K (sum_k L_k L_k)
        for i in range(d):
            for k in range(d):
                Kn -= 0.25 * s[i, i][:, None] * s[k, k][None, :] * K[np.ix_(r[i, i], r[k, k])]
        K = Kn
    return K


def sym_basis_full(n, d=3):
    """Columns: orthonormal symmetric basis vectors e_alpha in (R^d)^{(x)n} (d^n x N)."""
    import itertools
    al = mindices(n, d)
    idx = index_of(n, d)
    B = np.zeros((d ** n, len(al)))
    for w in itertools.product(range(d), repeat=n):
        cnt = [0] * d
        for c in w:
            cnt[c] += 1
        k = idx[tuple(cnt)]
        flat = 0
        for c in w:
            flat = flat * d + c
        B[flat, k] = 1.0
    B /= np.linalg.norm(B, axis=0)
    return B


def apply_Pn(X, n, a, d=3):
    """P_a^{(x)n}(X) for X a d^n x d^n matrix (real or complex)."""
    T = X.reshape([d] * (2 * n))
    for k in range(n):
        # move copy k's (row, col) indices to front
        T = np.moveaxis(T, [k, n + k], [0, 1])
        tr = np.trace(T, axis1=0, axis2=1)
        new = -0.25 * (T + np.swapaxes(T, 0, 1))
        for i in range(d):
            new[i, i] += a * tr
        T = np.moveaxis(new, [0, 1], [k, n + k])
    return T.reshape(d ** n, d ** n)


def V_full(x, y, n, a, d=3):
    X = np.outer(x, x.conj())
    M = apply_Pn(X, n, a, d)
    return np.real(np.vdot(y, M @ y))


def a_n(n):
    return (1 + 2 ** (1 - 1 / n)) / 4


def witness_sym(n, d=3):
    """Coordinates (in the monomial basis) of the witness x and y (symmetric) of Theorem 7.1."""
    B = sym_basis_full(n, d)
    p = np.zeros(d, complex); p[0] = 1 / sqrt(2); p[1] = 1j / sqrt(2)
    pn = p
    for _ in range(n - 1):
        pn = np.kron(pn, p)
    x = np.real((pn + pn.conj()) / sqrt(2))
    if n % 2 == 1:
        y = x.copy()
    else:
        y = np.real(1j * (pn - pn.conj()) / sqrt(2))
    return B.T @ x, B.T @ y


def Mx(K, x):
    N = len(x)
    K4 = K.reshape(N, N, N, N)
    return np.einsum('a,abcd,c->bd', x, K4, x)


def My(K, y):
    N = len(y)
    K4 = K.reshape(N, N, N, N)
    return np.einsum('b,abcd,d->ac', y, K4, y)


def alt_min(K, N, starts=50, iters=300, rng=None, x0=None):
    rng = np.random.default_rng(rng)
    best = (np.inf, None, None)
    for s in range(starts):
        if x0 is not None and s == 0:
            x = x0.copy()
        else:
            x = rng.standard_normal(N)
        x /= np.linalg.norm(x)
        val_old = np.inf
        for it in range(iters):
            w, v = np.linalg.eigh(Mx(K, x)); y = v[:, 0]
            w2, v2 = np.linalg.eigh(My(K, y)); x = v2[:, 0]
            val = w2[0]
            if abs(val_old - val) < 1e-15:
                break
            val_old = val
        if val < best[0]:
            best = (val, x.copy(), y.copy())
    return best
