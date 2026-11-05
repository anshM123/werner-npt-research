"""Exploration (floating point): on which quantities do G_z e_rho, G_z (0,1,...,1) and G_z v depend?
Fix s_t, s_1 (= s_h) and S = s_U, move the other U-weights, and watch coordinate h = 1."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cf import TypeI

np.set_printoptions(linewidth=200, precision=10, suppress=True)


def Gz(st, sU):
    w = np.sqrt(np.array([st] + list(sU), dtype=float))
    ti = TypeI(w, [0], F=float)
    M = np.array(ti.build())
    U = ti.U
    el = [ti.idx[('P', r, h)] for r in U for h in U] + [ti.idx[('Q',)]]
    tb = [ti.idx[('R', 0)]] + [ti.idx[('Z', 0, h)] for h in U]
    C = M[np.ix_(el, el)]
    B = M[np.ix_(el, tb)]
    A = M[np.ix_(tb, tb)]
    G = A - B.T @ np.linalg.solve(C, B)
    u = len(sU)
    Dg = np.diag([1.0] + [w[0] * w[1 + h] for h in range(u)])
    return Dg @ G @ Dg


if __name__ == "__main__":
    st, s1, S = 0.7, 0.4, 2.0
    for u in (3, 4):
        print(f'u = {u}')
        rng = np.random.default_rng(5)
        for it in range(4):
            rest = rng.uniform(0.1, 1, u - 1)
            rest = rest / rest.sum() * (S - s1)
            sU = [s1] + list(rest)
            G = Gz(st, sU)
            e = np.zeros(u + 1); e[0] = 1
            o = np.zeros(u + 1); o[1:] = 1
            print('  rest', rest, '\n    G e_rho =', G @ e, '\n    G (0,1..1) =', G @ o, '\n    G_00 =', G[0, 0])
