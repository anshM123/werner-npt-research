"""Three-replica SDP relaxation: min q_n(Pi) over rank-two projectors whose range is t-uniform.
Uses: T = twirl(Pi^{(x)3}) >= 0 (copy-symmetrised blocks X_j >= 0), global replica symmetry, Lambda^3(L) = 0
(T P_sign = 0), global moments Tr[T tau] = 2^{cyc(tau)}, A_j = B_j = 0 for 1 <= j <= t, B_j >= A_j."""
import sys
import numpy as np
import cvxpy as cp
from rep3 import (PERMS, E, T12, C1, ncycles, blocks, product_functional, sign_projector, group_action)


def build(n, t, shadow=True):
    blks = blocks(n)
    X = {bk: cp.Variable((bk[3] + 1, bk[3] + 1), symmetric=True) for bk in blks}
    cons = [X[bk] >> 0 for bk in blks]

    def lin(funcs, k):
        terms = []
        for bk, mats in funcs.items():
            if k < len(mats) and np.abs(mats[k]).max() > 0:
                terms.append(cp.sum(cp.multiply(X[bk], mats[k])))
        return cp.sum(terms) if terms else 0

    for tau in PERMS:
        f = product_functional(n, {tau: np.array([1.0])})
        cons.append(lin(f, 0) == 2.0 ** ncycles(tau))
    for bk in blks:
        a, m, b, j2 = bk
        P = sign_projector(m, j2, b)
        if np.abs(P).max() > 1e-12:
            cons.append(cp.trace(X[bk] @ P) == 0)
        for g in (T12, C1):
            G = group_action(m, j2, g)
            cons.append(X[bk] @ G == G @ X[bk])
    fA = product_functional(n, {E: np.array([1.0, -1.0]), T12: np.array([0.0, 3.0])})
    fB = product_functional(n, {E: np.array([0.0, 3.0]), T12: np.array([1.0, -1.0])})
    Aj = [lin(fA, k) / 8 for k in range(n + 1)]
    Bj = [lin(fB, k) / 4 for k in range(n + 1)]
    for j in range(1, t + 1):
        cons += [Aj[j] == 0, Bj[j] == 0]
    for j in range(2, n + 1):
        cons.append(Bj[j] - Aj[j] >= 0)
        cons.append(Aj[j] >= 0)
    if shadow:
        fS = product_functional(n, {E: np.array([0.5, 0.5]), T12: np.array([0.5, -0.5])})
        for k in range(n + 1):
            cons.append(lin(fS, k) >= 0)
    fO = product_functional(n, {T12: np.array([2.0]), E: np.array([-1.0])})
    obj = lin(fO, 0) / 2 / 2 ** n
    return cp.Problem(cp.Minimize(obj), cons), Aj, Bj


if __name__ == '__main__':
    n0, n1 = int(sys.argv[1]), int(sys.argv[2])
    t = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    for n in range(n0, n1 + 1):
        prob, Aj, Bj = build(n, t)
        try:
            val = prob.solve(solver=cp.CLARABEL)
            print(f'n={n:2d} t={t}: SDP min q = {val:+.6f}  status {prob.status}', flush=True)
        except Exception as e:
            print(f'n={n:2d} t={t}: solver error {e}', flush=True)
