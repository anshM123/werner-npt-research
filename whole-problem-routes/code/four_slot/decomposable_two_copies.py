"""Decisive necessary test for EVERY frame-based tensor-stable certificate (sandwich chains, GM trees, any DAG of
partial-transpose block-PSD conditions on per-site operators):  such a certificate at parameter t implies
Tr[P (A^{(x)n} - B^{(x)n})] >= 0 for every 4-party state P that is PPT in every frame, i.e. A^{(x)n} - B^{(x)n}
is 4-party decomposable for every n.  Here we test n = 1 and n = 2 by SDP, exactly reduced by the U(3)^{(x)4}
symmetry at each site: frame-S operators live in A_S = span{P_pi^{G_S}} (23-dim *-algebra); PSD is tested in the
irreducible blocks of its regular representation (faithful *-rep for the trace inner product).
mu* = max{mu : A^{(x)n}-B^{(x)n} = sum_S X_S, X_S^{G_S} >= mu I}.  mu* >= 0 (to tolerance) <=> decomposable.
"""
import sys, time
import numpy as np
import cvxpy as cp
from common import AB, ptrans, subsets, all_perms

perms, Pm = all_perms()
NP = len(perms)
FR = [S for S in subsets() if 0 not in S]  # 8 frames up to complement
rng = np.random.default_rng(7)


def onb(mats):
    V = np.array([m.ravel() for m in mats])
    U, s, Wt = np.linalg.svd(V, full_matrices=False)
    k = (s > 1e-9 * s[0]).sum()
    return [Wt[i].reshape(81, 81) for i in range(k)]


def irreps(S):
    E = [ptrans(P, S) for P in Pm]
    F = onb(E)
    k = len(F)
    Fv = np.array([f.ravel() for f in F])
    L = lambda X: np.array([[np.sum(F[a] * (X @ F[b])) for b in range(k)] for a in range(k)])
    R = lambda X: np.array([[np.sum(F[a] * (F[b] @ X)) for b in range(k)] for a in range(k)])
    Z = sum(rng.normal() * (e + e.T) / 2 for e in E)
    w, V = np.linalg.eigh(R(Z))
    groups, cur = [], [0]
    for i in range(1, k):
        if abs(w[i] - w[cur[-1]]) < 1e-7:
            cur.append(i)
        else:
            groups.append(cur); cur = [i]
    groups.append(cur)
    Qs = [V[:, g] for g in groups]
    LE = [L(e) for e in E]
    reps = [[Q.T @ Le @ Q for Le in LE] for Q in Qs]   # reps[i][pi] (m_i x m_i)
    # check: images of A-B style elements are symmetric when they should be
    return reps, F, E


def setup():
    data = {}
    for S in FR:
        reps, F, E = irreps(S)
        data[S] = reps
        print(f"frame {S}: block sizes {[r[0].shape[0] for r in reps]}", flush=True)
    F0 = onb(Pm)
    T = np.array([[np.sum(f * P) for P in Pm] for f in F0])  # coefficient map perms -> ONB of frame ()
    return data, F0, T


def solve(t, n, data, F0, T):
    A, B = AB(t)
    a = np.array([np.sum(f * A) for f in F0]); b = np.array([np.sum(f * B) for f in F0])
    mu = cp.Variable()
    cons = []
    if n == 1:
        target = a - b
        ys = {S: cp.Variable(NP) for S in FR}
        cons.append(sum(T @ ys[S] for S in FR) == target)
        for S in FR:
            for rep in data[S]:
                m = rep[0].shape[0]
                K = np.array([r.ravel(order='F') for r in rep]).T  # (m*m) x NP
                M = cp.reshape(K @ ys[S], (m, m), order='F')
                cons.append((M + M.T) / 2 - mu * np.eye(m) >> 0)
    else:
        target = np.outer(a, a) - np.outer(b, b)
        ys = {S: cp.Variable((NP, NP)) for S in FR}
        cons.append(sum(T @ ys[S] @ T.T for S in FR) == target)
        for S in FR:
            reps = data[S]
            for ri in reps:
                for rj in reps:
                    mi, mj = ri[0].shape[0], rj[0].shape[0]
                    m = mi * mj
                    K = np.array([np.kron(ri[p], rj[q]).ravel(order='F') for p in range(NP) for q in range(NP)]).T
                    M = cp.reshape(K @ cp.vec(ys[S].T, order='F'), (m, m), order='F')  # row-major vec of y
                    cons.append((M + M.T) / 2 - mu * np.eye(m) >> 0)
    for S in FR:
        cons.append(cp.norm(cp.vec(ys[S], order='F'), 'inf') <= 1e3)
    prob = cp.Problem(cp.Maximize(mu), cons)
    for solver in (cp.CLARABEL, cp.SCS):
        try:
            prob.solve(solver=solver) if solver == cp.CLARABEL else prob.solve(solver=solver, eps=1e-9, max_iters=100000)
            if prob.status in ('optimal', 'optimal_inaccurate'):
                return mu.value, prob.status + ('' if solver == cp.CLARABEL else '[SCS]')
        except Exception as e:
            print('  solver error', solver, e, flush=True)
    return None, prob.status


if __name__ == '__main__':
    data, F0, T = setup()
    for t in [0.5, 0.55]:
        t0 = time.time(); mu, st = solve(t, 1, data, F0, T)
        print(f"n=1 t={t}: mu*={mu} ({st}) [{time.time()-t0:.1f}s]", flush=True)
    for t in [1/3, 0.34, 0.36, 0.40, 0.45, 0.5]:
        t0 = time.time(); mu, st = solve(t, 2, data, F0, T)
        print(f"n=2 t={t:.4f}: mu*={mu} ({st}) [{time.time()-t0:.1f}s]", flush=True)
