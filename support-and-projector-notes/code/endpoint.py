"""Shared helper: the endpoint form q_n(X, Y) = <X, Phi_n(Y)>, Phi_n = (id - t Tr(.) I)^{(x) n}, local dimension 3."""
import numpy as np

d = 3


def Phi(Y, n, t=0.5):
    T = Y.reshape([d] * (2 * n))
    for s in range(n):
        tr = np.trace(T, axis1=s, axis2=n + s)
        full = np.moveaxis(np.einsum('ab,...->...ab', np.eye(d), tr), [2 * n - 2, 2 * n - 1], [s, n + s])
        T = T - t * full
    return T.reshape(d ** n, d ** n)


def q(X, n, Y=None, t=0.5):
    """Sesquilinear in (X, Y), conjugate-linear in X; q(X, n) is the quadratic form."""
    return np.vdot(X, Phi(X if Y is None else Y, n, t))
