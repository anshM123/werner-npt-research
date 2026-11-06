"""Floating-point helpers: the Schur complement G on the letter-2 states (rho_t, zeta_{t,h}) of a type I block of
(R_inf)_{k-1}, built from the closed form cf.TypeI (units of 2Y), and its frame (rho_t, z_t, zeta'_t), where
zeta_t = z_t w_U/|w_U| + zeta'_t with zeta'_t orthogonal to w_U.  Exploration only."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cf import TypeI


def Gmat(sT, sU):
    """sT, sU: squared weights of T and U.  Copies ordered T first, then U.  Returns G (ordered per t: rho_t, then
    zeta_{t,h}, h in U) and the full block."""
    w = np.sqrt(np.array(list(sT) + list(sU), dtype=float))
    m = len(sT)
    T = list(range(m))
    ti = TypeI(w, T, F=float)
    M = np.array(ti.build())
    n = ti.n
    U = ti.U
    r0 = U[0]
    tb = []
    for t in T:
        tb.append(ti.idx[('R', t)])
        for h in U:
            tb.append(ti.idx[('Z', t, h)])
    el = [i for i in range(n) if i not in tb and not (ti.labels[i][0] == 'P' and ti.labels[i][1] == r0)]
    A11 = M[np.ix_(el, el)]
    A12 = M[np.ix_(el, tb)]
    A22 = M[np.ix_(tb, tb)]
    G = A22 - A12.T @ np.linalg.solve(A11, A12)
    return (G + G.T) / 2


def frame(sT, sU, rescale=True):
    """Block-diagonal change of coordinates P with zeta_t = c z_t w_U/|w_U| + B zeta'_t (B orthonormal basis of
    w_U-perp), c = sqrt(1 + S/sigma) if rescale.  Returns P (so that Ghat = P^T G P)."""
    u = len(sU)
    m = len(sT)
    d = u + 1
    wU = np.sqrt(np.array(sU, dtype=float))
    S = sum(sU)
    sig = sum(sT)
    wh = wU / np.linalg.norm(wU)
    c = np.sqrt(1 + S / sig) if rescale else 1.0
    Pt = np.zeros((d, d))
    Pt[0, 0] = 1
    Pt[1:, 1] = c * wh
    if u > 1:
        Bm = np.linalg.svd(wh[None, :])[2][1:].T
        Pt[1:, 2:] = Bm
    P = np.zeros((m * d, m * d))
    for t in range(m):
        P[t * d:(t + 1) * d, t * d:(t + 1) * d] = Pt
    return P


def Ghat(sT, sU, rescale=True):
    G = Gmat(sT, sU)
    P = frame(sT, sU, rescale)
    return P.T @ G @ P
