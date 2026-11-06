"""Exact (rational) evaluation of G_z for m = 1 at rational weights w (so s = w^2 is rational).

The type I block of (R_inf)_{k-1} with T = {t} = {0} and U = {1, ..., u} is built by cf.TypeI (units of 2Y) with
python-flint fmpq entries.  Gauge: the coordinates x_t = (P_{t,h})_h are deleted.  This is a complement of the kernel
span{p (x) f_h} (the kernel vector of c has P_{t,h}-coordinate w_t c_h), so the Schur complement onto the letter-2 states
does not depend on this choice (it is also compared below with the gauge P_{r0,h} deleted, r0 = 1).  The eliminated
block is then the (C_u) block on (X_U, q), positive definite by Theorem C.

G_z = Dg G Dg with Dg = diag(1, w_t w_h): the coordinates rho and z_h with zeta_h = w_t w_h z_h.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from cf import TypeI

Q = flint.fmpq


def full_block(wt, wU):
    w = [Q(wt)] + [Q(x) for x in wU]
    ti = TypeI(w, [0], F=Q)
    M = ti.build()
    return w, ti, M


def schur_G(w, ti, M, gauge='t'):
    U = ti.U
    tb = [ti.idx[('R', 0)]] + [ti.idx[('Z', 0, h)] for h in U]
    if gauge == 't':
        el = [ti.idx[('P', r, h)] for r in U for h in U] + [ti.idx[('Q',)]]
    else:
        r0 = U[0]
        el = [ti.idx[('P', r, h)] for r in [0] + U if r != r0 for h in U] + [ti.idx[('Q',)]]
    C = flint.fmpq_mat([[M[i][j] for j in el] for i in el])
    B = flint.fmpq_mat([[M[i][j] for j in tb] for i in el])
    A = flint.fmpq_mat([[M[i][j] for j in tb] for i in tb])
    Y = C.solve(B)
    G = A - B.transpose() * Y
    return G, C


def Gz_exact(wt, wU, gauge='t', return_C=False):
    w, ti, M = full_block(wt, wU)
    G, C = schur_G(w, ti, M, gauge)
    u = len(wU)
    d = [Q(1)] + [w[0] * w[1 + h] for h in range(u)]
    Gz = flint.fmpq_mat([[G[i, j] * d[i] * d[j] for j in range(u + 1)] for i in range(u + 1)])
    if return_C:
        return Gz, C
    return Gz


def mvec(st, S, u, c1=4, c2=3):
    """v = (4 s_t (s_t + S), 3 (2 s_t + S), ..., 3 (2 s_t + S))."""
    return [c1 * st * (st + S)] + [c2 * (2 * st + S)] * u


def ldl_pivots(Mq):
    n = Mq.nrows()
    A = [[Mq[i, j] for j in range(n)] for i in range(n)]
    piv = []
    for p in range(n):
        pv = A[p][p]
        piv.append(pv)
        if pv == 0:
            return piv
        for i in range(p + 1, n):
            f = A[i][p] / pv
            if f != 0:
                for j in range(p + 1, n):
                    A[i][j] -= f * A[p][j]
    return piv
