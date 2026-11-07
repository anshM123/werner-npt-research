"""Blocks of the compressed words on Sym^n (x) Sym^n, spin by spin, in the orthonormal coupled basis of cpl.py.

Symbols (coherent states u^n (x) w^n, u'^n (x) w'^n):  p = (u.u')(w.w'),  q = (u.w')(w.u'),  r = (u.w)(u'.w').
M(k,s,t) := compression of ONE arrangement of the word I^k F^s Omega^t; symbol p^k q^s r^t.
Factorisation (Fock picture, m = n - t):
    M(k,s,t) = (m!/n!)^2 (a+.b+)^t D^{(m)}_s (a.b)^t,   D^{(m)}_s = sum_j omega_j^{(m)}(s) P_j^{(m)},
P_j^{(m)} = projector onto the GL(3)-type (2m-j, j) part of Sym^m (x) Sym^m = eigenspace of E_ab E_ba with eigenvalue
(m-j)(m-j+1).  Hence M(k,s,t) = sum_j omega_j(s) Zt_{t,j} with the positive semidefinite
    Zt_{t,j} = (m!/n!)^2 (a+.b+)^t P_j^{(m)} (a.b)^t = Y_{t,j}^T Y_{t,j},   Y_{t,j} = rows of V_m^T Xt_t with label j.
"""
import math
import numpy as np
from cpl import basis, ab, Eba, omega, mu


def level_data(n, L):
    """Returns list over t = 0..n of (Y_t, labels_t): Y_t (rows: eigenvectors of E_ab E_ba at level m = n - t,
    columns: (n,n) spin-L basis), labels_t[i] = j."""
    d = len(basis(n, n, L))
    X = np.eye(d)
    out = []
    for t in range(0, n + 1):
        m = n - t
        if t > 0:
            X = (ab(m + 1, L) / (m + 1)) @ X
        dm = len(basis(m, m, L))
        if dm == 0:
            out.append((np.zeros((0, d)), np.zeros(0, int)))
            continue
        E = Eba(m, L)
        C = E.T @ E
        ev, V = np.linalg.eigh(C)
        allowed = np.array([(m - j) * (m - j + 1) for j in range(m + 1)], float)
        labels = np.array([int(np.argmin(np.abs(allowed - e))) for e in ev])
        assert np.abs(allowed[labels] - ev).max() < 1e-6 * max(1.0, allowed.max())
        out.append((V.T @ X, labels))
    return out


def word_block(ld, n, k, s, t):
    """Block of M(k,s,t) from level data ld = level_data(n, L)."""
    m = n - t
    assert k + s == m
    Y, lab = ld[t]
    w = np.array([float(omega(m, j, s)) for j in lab])
    return (Y.T * w) @ Y


def Zt_block(ld, t, j):
    Y, lab = ld[t]
    R = Y[lab == j]
    return R.T @ R
