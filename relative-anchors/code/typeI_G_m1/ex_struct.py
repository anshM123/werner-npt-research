"""Exploration (floating point, no claim rests on it): m = 1, the minimiser X* of the eliminated block (gauge x_t = 0)
for the letter-2 data v = (rho, zeta) = (1, beta w_t w_U).  Prints the structure of X*_U in scaled coordinates."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cf import TypeI

np.set_printoptions(linewidth=200, precision=6, suppress=True)


def blocks(st, sU):
    w = np.sqrt(np.array([st] + list(sU), dtype=float))
    ti = TypeI(w, [0], F=float)
    M = np.array(ti.build())
    U = ti.U
    el = [ti.idx[('P', r, h)] for r in U for h in U] + [ti.idx[('Q',)]]
    tb = [ti.idx[('R', 0)]] + [ti.idx[('Z', 0, h)] for h in U]
    C = M[np.ix_(el, el)]
    B = M[np.ix_(el, tb)]
    A = M[np.ix_(tb, tb)]
    return w, ti, M, C, B, A, el, tb


rng = np.random.default_rng(1)
if __name__ == "__main__":
  for u in (2, 3, 4):
    st = float(np.exp(rng.uniform(-1, 1)))
    sU = list(np.exp(rng.uniform(-1.5, 1.5, u)))
    w, ti, M, C, B, A, el, tb = blocks(st, sU)
    S = sum(sU)
    beta = 3 * (2 * st + S) / (4 * st * (st + S))
    vact = np.array([1.0] + [beta * w[0] * w[1 + h] for h in range(u)])
    xs = -np.linalg.solve(C, B @ vact)
    X = xs[:-1].reshape(u, u)
    q = xs[-1]
    wU = w[1:]
    Xs = X / np.outer(wU, wU)
    print(f'u = {u}, s_t = {st:.4f}, s_U = {np.array(sU)}')
    print('X* (scaled x_rh = X_rh/(w_r w_h)):')
    print(Xs)
    print('symmetric?', np.max(np.abs(Xs - Xs.T)))
    print('diag of X (unscaled):', np.diag(X))
    off = Xs.copy()
    np.fill_diagonal(off, np.nan)
    print('q =', q)
    # compare G v with closed G
    G = A - B.T @ np.linalg.solve(C, B)
    Dg = np.diag([1.0] + [w[0] * w[1 + h] for h in range(u)])
    Gz = Dg @ G @ Dg
    vz = np.array([1.0] + [beta] * u)
    print('G_z v =', Gz @ vz)
    print()
