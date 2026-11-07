"""Restricted margin SDP at tau' = 1: only the pairs (s,t), 1 <= t < s, with k = n - s - t <= K0 are free; all other pairs
keep the greedy value y = 1 (the pairs (s,0) are forced to 1 by the kernel Sym^{2n}).  Maximise the absolute least
eigenvalue lam of Q on the complement of Sym^{2n} (unscaled).  lam > 0 means the decomposition bound is exactly 4^{-n}.
Blocks are split by the input/output swap S (all words commute with S).
usage: python restricted_sdp.py K0 ybound n1 n2 ..."""
import sys, time, pickle
import numpy as np, cvxpy as cp
from scipy.linalg import null_space
from wsdp import Blocks, Q0_poly, F_poly, pairs, add
from cpl import swap_sign

def sym_split(n, L):
    S = swap_sign(n, L)
    ev, V = np.linalg.eigh((S + S.T) / 2)
    return [V[:, ev > 0], V[:, ev < 0]]

def setup(n, K0):
    B = Blocks(n)
    prs = pairs(n)
    free = [p for p in prs if p[1] >= 1 and n - p[0] - p[1] <= K0]
    P0 = Q0_poly(n)
    for p in prs:
        if p not in free:
            P0 = add(P0, F_poly(n, *p))
    Fpolys = [F_poly(n, *p) for p in free]
    data = []
    for L in B.Ls:
        Q0 = B.block(P0, L); Fs = [B.block(Fp, L) for Fp in Fpolys]
        d = Q0.shape[0]
        if L % 2 == 0:
            Y0, lab0 = B.ld[L][0]
            v = Y0[lab0 == 0][0]; v = v / np.linalg.norm(v)
            assert np.linalg.norm(Q0 @ v) < 1e-6 * max(1, np.abs(Q0).max()), 'kernel'
            Pv = np.eye(d) - np.outer(v, v)
        else:
            Pv = np.eye(d)
        for U in sym_split(n, L):
            U = U - np.outer(np.zeros(d), np.zeros(U.shape[1]))
            W = Pv @ U
            # orthonormal basis of the range of W
            if W.shape[1] == 0:
                continue
            uu, sv, _ = np.linalg.svd(W, full_matrices=False)
            W = uu[:, sv > 1e-9]
            if W.shape[1] == 0:
                continue
            data.append((L, W.T @ Q0 @ W, [W.T @ F @ W for F in Fs]))
    return free, data

def solve(n, K0, ybound=20.0):
    t0 = time.time()
    free, data = setup(n, K0)
    t1 = time.time()
    y = cp.Variable(len(free)); lam = cp.Variable()
    cons = [cp.abs(y) <= ybound]
    for (L, Q0, Fs) in data:
        M = Q0 + (sum(y[i] * Fs[i] for i in range(len(free))) if free else 0)
        cons.append((M + M.T) / 2 - lam * np.eye(Q0.shape[0]) >> 0)
    prob = cp.Problem(cp.Maximize(lam), cons)
    try:
        prob.solve(solver='CLARABEL')
    except Exception as e:
        prob.solve(solver='SCS', eps=1e-9, max_iters=300000)
    yv = y.value
    # independent recheck of the least eigenvalue at yv
    lmin = np.inf; Lat = None
    for (L, Q0, Fs) in data:
        M = Q0 + sum(yv[i] * Fs[i] for i in range(len(free)))
        e = np.linalg.eigvalsh((M + M.T) / 2)[0]
        if e < lmin:
            lmin, Lat = e, L
    return dict(status=prob.status, lam=lam.value, lmin=lmin, Lat=Lat, y=yv, free=free, tsetup=t1 - t0, tsolve=time.time() - t1)

if __name__ == '__main__':
    K0 = int(sys.argv[1]); yb = float(sys.argv[2])
    for n in [int(x) for x in sys.argv[3:]]:
        r = solve(n, K0, yb)
        T = dict(zip(r['free'], r['y']))
        print(f"n={n} K0={K0}: {r['status']}, margin {r['lam']:.5f}, recomputed least eigenvalue {r['lmin']:.5f} (L={r['Lat']}); "
              f"{len(r['free'])} free pairs; setup {r['tsetup']:.1f}s solve {r['tsolve']:.1f}s", flush=True)
        for k in range(K0 + 1):
            row = [(p, T[p]) for p in r['free'] if n - p[0] - p[1] == k]
            print(f"   k={k}: " + ' '.join(f"{p[0]}-{p[1]}:{v:.3f}" for p, v in row), flush=True)
        pickle.dump(dict(n=n, K0=K0, **r), open(f'restr_K{K0}_n{n}.pkl', 'wb'))
