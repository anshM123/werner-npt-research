"""Two-replica (degree-four) relaxation of q_N(t, P) on real rank-two orthogonal projectors P = vv^T + ww^T whose range
is balanced (every one-copy compression equals (Tr H/d) I_2), symmetrised over permutations of the copies.
Variables: M_k[(a,b),(c,e)] = <Tr_K(a b^T), Tr_K(c e^T)> for |K| = k and a, b, c, e in {v, w}.
Constraints, each valid for real balanced projectors: M_k PSD; M_0 = delta_ac delta_be; balance M_{N-1} = delta_ab delta_ce / d;
the crossing M_{N-k}[(a,b),(c,e)] = M_k[(e,b),(c,a)]; M[(p,q),(r,s)] = M[(s,r),(q,p)]; M[(a,b),(c,e)] = M[(b,a),(e,c)];
the Gram matrices of the Werner sectors Pi_T(a (x) b) and of the identity components of a b^T are PSD.
Prints the minimum of q_N(t, P)/||P||_F^2 over the relaxation.  At t = 1/2 the values coincide with the weight-enumerator
bounds of Proposition 19 (3/4, 5/16, 1239/10688 per unit ||P||_F^2 at N = 3, 5, 7, and the negative value at N = 9).
usage: python two_replica_balanced.py t d N_1 N_2 ..."""
import sys
from math import comb
import numpy as np
import cvxpy as cp

t = float(sys.argv[1]); d = int(sys.argv[2]); Ns = [int(x) for x in sys.argv[3:]]
pairs = [(a, b) for a in range(2) for b in range(2)]
idx = {p: i for i, p in enumerate(pairs)}


def pm(fun):
    P = np.zeros((16, 16))
    for (a, b) in pairs:
        for (c, e) in pairs:
            (p, q), (r, s) = fun(a, b, c, e)
            P[idx[(c, e)] * 4 + idx[(a, b)], idx[(r, s)] * 4 + idx[(p, q)]] = 1      # column-major vec
    return P


P_cross = pm(lambda a, b, c, e: ((e, b), (c, a)))
P_adj = pm(lambda a, b, c, e: ((e, c), (b, a)))
P_tr = pm(lambda a, b, c, e: ((b, a), (e, c)))
P_F = pm(lambda a, b, c, e: ((a, c), (e, b)))        # <a(x)b|F_L|c(x)e> = M_{L^c}[(a,c),(e,b)]


def re(M, P):
    if isinstance(M, np.ndarray):
        return (P @ M.reshape(-1, order='F')).reshape((4, 4), order='F')
    return cp.reshape(P @ cp.vec(M, order='F'), (4, 4), order='F')


def solve(N):
    M0 = np.array([[float(a == c and b == e) for (c, e) in pairs] for (a, b) in pairs])
    Mbal = np.array([[float(a == b and c == e) / d for (c, e) in pairs] for (a, b) in pairs])
    M = {0: M0, N: re(M0, P_cross), N - 1: Mbal, 1: re(Mbal, P_cross)}
    cons = []
    for k in range(2, N - 1):
        if k in M:
            continue
        V = cp.Variable((4, 4), symmetric=True)
        M[k] = V
        if N - k != k:
            M[N - k] = re(V, P_cross)
        else:
            cons.append(V == re(V, P_cross))
        cons += [V == re(V, P_adj), V == re(V, P_tr)]
    for k in range(N + 1):
        if not isinstance(M[k], np.ndarray):
            cons.append((M[k] + M[k].T) / 2 >> 0)
    for j in range(N + 1):
        S = 0
        for l in range(N + 1):
            coef = sum(comb(j, i) * comb(N - j, l - i) * (-1) ** i for i in range(max(0, l - (N - j)), min(j, l) + 1))
            if coef:
                S = S + coef * re(M[N - l], P_F)
        S = S / 2 ** N
        if isinstance(S, np.ndarray):
            assert np.linalg.eigvalsh((S + S.T) / 2).min() > -1e-9
        else:
            cons.append((S + S.T) / 2 >> 0)
    for k in range(N + 1):
        Ck = sum(comb(N - k, l - k) * (-1) ** (l - k) * d ** (-l) * M[l] for l in range(k, N + 1))
        if not isinstance(Ck, np.ndarray):
            cons.append((Ck + Ck.T) / 2 >> 0)
    diag = [idx[(0, 0)], idx[(1, 1)]]
    q = sum(comb(N, k) * (-t) ** k * sum(M[k][i, j] for i in diag for j in diag) for k in range(N + 1))
    if not cons:
        return float(q) / 2, 'fixed'
    prob = cp.Problem(cp.Minimize(q), cons)
    prob.solve(solver='CLARABEL')
    return prob.value / 2, prob.status


for N in Ns:
    val, st = solve(N)
    print(f"t={t} d={d} N={N}: relaxation min q_N(t,P)/||P||_F^2 = {val:+.6e}  [{st}]", flush=True)
