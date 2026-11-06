"""Floating-point helpers (exploration only): the type I block of (R_inf)_{k-1} with T = {0} (m = 1), U = {1..u},
in the gauge x_t = 0 (states P_{0,h} deleted), from the closed form cf.TypeI (units of 2Y)."""
import numpy as np
from cf import TypeI


def block(st, sU, scaled=True):
    """Returns (Phi, labels): the full block on (X_U (row-major r,h), q, rho, zeta_h) in gauge x_t = 0.
    scaled: X_rh = w_r w_h x_rh, zeta_h = w_t w_h z_h."""
    w = np.sqrt(np.array([st] + list(sU), dtype=float))
    ti = TypeI(w, [0], F=float)
    M = np.array(ti.build())
    U = ti.U
    el = [ti.idx[('P', r, h)] for r in U for h in U] + [ti.idx[('Q',)]]
    tb = [ti.idx[('R', 0)]] + [ti.idx[('Z', 0, h)] for h in U]
    order = el + tb
    P = M[np.ix_(order, order)]
    if scaled:
        sc = np.array([w[ti.labels[i][1]] * w[ti.labels[i][2]] if ti.labels[i][0] in ('P', 'Z') else 1.0
                       for i in order])
        P = P * np.outer(sc, sc)
    return P, [ti.labels[i] for i in order], len(el)


def Gz(st, sU):
    P, lab, ne = block(st, sU)
    C = P[:ne, :ne]; B = P[:ne, ne:]; A = P[ne:, ne:]
    G = A - B.T @ np.linalg.solve(C, B)
    return (G + G.T) / 2, C, B, A


def vvec(st, S, u):
    return np.array([4 * st * (st + S)] + [3 * (2 * st + S)] * u)
