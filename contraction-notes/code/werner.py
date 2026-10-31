"""Shared helpers. Operators on (C^d)^{(x)k} are arrays of shape (d^k, d^k), rows then columns.

    q_k(t, X) = sum_J (-t)^{|J|} ||Tr_J X||_2^2 = <X, B(X)>,   B = prod_s (1 - t E_s),   E_s X = I_s (x) Tr_s X.

Crossed pairing: for rank-one operators q_k(t; |x><y|, |u><v|) = <x (x) v| W |u (x) y> with W = prod_s (I - t F_s),
F_s exchanging the s-th factors of two copies."""
import itertools
import numpy as np

d = 3


def applyB(X, k, t):
    T = X.reshape([d] * (2 * k))
    for s in range(k):
        tr = np.trace(T, axis1=s, axis2=k + s)
        E = np.multiply.outer(tr, np.eye(d))
        E = np.moveaxis(E, [2 * k - 2, 2 * k - 1], [s, k + s])
        T = T - t * E
    return T.reshape(d ** k, d ** k)


def q(X, k, t):
    return np.vdot(X, applyB(X, k, t)).real


def q_pair(X, Y, k, t):
    """sesquilinear form q_k(t; X, Y), conjugate-linear in X"""
    return np.vdot(X, applyB(Y, k, t))


def q_direct(X, k, t):
    """the defining sum over subsets J, for checking applyB"""
    T = X.reshape([d] * (2 * k)); tot = 0.0
    for J in itertools.product([0, 1], repeat=k):
        A, mm = T, k
        for i in reversed(range(k)):
            if J[i]:
                A = np.trace(A, axis1=i, axis2=i + mm); mm -= 1
        tot += (-t) ** sum(J) * np.vdot(A, A).real
    return tot


def W_apply(v, k, t):
    """W = prod_s (I - t F_s) applied to v in (C^d)^{(x)k} (x) (C^d)^{(x)k}"""
    T = v.reshape([d] * (2 * k))
    for s in range(k):
        T = T - t * np.swapaxes(T, s, k + s)
    return T.reshape(-1)


def contract(C, k, s, w, wp):
    """<w|_s C |w'>_s"""
    T = C.reshape([d] * (2 * k))
    T = np.tensordot(np.conj(w), T, axes=([0], [s]))
    T = np.tensordot(T, wp, axes=([k - 1 + s], [0]))
    return T.reshape(d ** (k - 1), d ** (k - 1))


def unpack(x, D):
    Y = (x[:2 * D] + 1j * x[2 * D:4 * D]).reshape(D, 2)
    Z = (x[4 * D:6 * D] + 1j * x[6 * D:]).reshape(D, 2)
    return Y, Z


def pack(Y, Z):
    return np.concatenate([Y.real.ravel(), Y.imag.ravel(), Z.real.ravel(), Z.imag.ravel()])
