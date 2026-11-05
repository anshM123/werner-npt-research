"""Exploration (floating point): is the off-diagonal part of the minimiser X* (scaled) for the source v separable
(x_rh = a_r + b_h, or of low rank)?  Gauge x_t = 0."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cf import TypeI

np.set_printoptions(linewidth=200, precision=8, suppress=True)


def solve(st, sU, vz):
    w = np.sqrt(np.array([st] + list(sU), dtype=float))
    ti = TypeI(w, [0], F=float)
    M = np.array(ti.build())
    U = ti.U
    u = len(U)
    el = [ti.idx[('P', r, h)] for r in U for h in U] + [ti.idx[('Q',)]]
    tb = [ti.idx[('R', 0)]] + [ti.idx[('Z', 0, h)] for h in U]
    C = M[np.ix_(el, el)]
    B = M[np.ix_(el, tb)]
    Dg = np.array([1.0] + [w[0] * w[1 + h] for h in range(u)])
    xs = -np.linalg.solve(C, B @ (Dg * vz))
    X = xs[:-1].reshape(u, u)
    wU = w[1:]
    return X / np.outer(wU, wU), xs[-1], X


if __name__ == "__main__":
    rng = np.random.default_rng(7)
    u = 5
    st = 0.9
    sU = list(np.exp(rng.uniform(-1.5, 1.5, u)))
    S = sum(sU)
    beta = 3 * (2 * st + S) / (4 * st * (st + S))
    for name, vz in [('v', np.array([1.0] + [beta] * u)), ('e_rho', np.array([1.0] + [0.0] * u)),
                     ('merged', np.array([0.0] + [1.0] * u)), ('e_z1', np.array([0.0, 1.0] + [0.0] * (u - 1)))]:
        Xs, q, X = solve(st, sU, vz)
        # additive separable fit on off-diagonal entries
        rows = []
        rhs = []
        for r in range(u):
            for h in range(u):
                if r != h:
                    a = np.zeros(2 * u)
                    a[r] = 1
                    a[u + h] = 1
                    rows.append(a)
                    rhs.append(Xs[r, h])
        A = np.array(rows)
        sol, res, rk, sv = np.linalg.lstsq(A, np.array(rhs), rcond=None)
        resid = np.array(rhs) - A @ sol
        off = Xs.copy()
        np.fill_diagonal(off, 0)
        print(f'{name}: additive-separable residual {np.abs(resid).max():.3e}; asym {np.abs(Xs - Xs.T).max():.3e}; '
              f'sing. values of off-diag scaled part {np.linalg.svd(off, compute_uv=False)}')
