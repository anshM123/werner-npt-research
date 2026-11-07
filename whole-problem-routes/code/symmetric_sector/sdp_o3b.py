import numpy as np, cvxpy as cp, time, sys, pickle
from bos import *
from o3red import O3Red

def setup(n):
    R = O3Red(n)
    G = R.gamma_matrix()
    return R, G

def dec_sdp(R, G, K, tfix=None, solver='CLARABEL'):
    labels = R.labels; pos = R.pos
    ut = [i for i, (L, r, s) in enumerate(labels) if r <= s]
    Qs = {L: cp.Variable((R.blocks[L]['m'],)*2, symmetric=True) for L in R.Lvals}
    t = cp.Variable() if tfix is None else tfix
    qvec = cp.hstack([Qs[L][r, s] for (L, r, s) in labels])
    kvec = np.array([R.block(K, L)[r, s] for (L, r, s) in labels])
    ivec = np.array([1.0 if r == s else 0.0 for (L, r, s) in labels])
    Gu = G[ut, :]
    cons = [kvec[ut] - t*ivec[ut] - qvec[ut] - Gu @ qvec == 0] + [Qs[L] >> 0 for L in R.Lvals]
    if tfix is None:
        prob = cp.Problem(cp.Maximize(t), cons)
    else:
        prob = cp.Problem(cp.Minimize(0), cons)
    try:
        prob.solve(solver=solver)
    except Exception as e:
        return None, 'fail '+str(e)[:60], None
    tv = t.value if tfix is None else tfix
    return tv, prob.status, {L: Qs[L].value for L in R.Lvals}

if __name__ == '__main__':
    for n in [int(s) for s in sys.argv[1:]]:
        t0 = time.time()
        R, G = setup(n)
        print(f"n={n}: N={R.N}, dim_inv={R.dim_inv}, m_L={[R.blocks[L]['m'] for L in R.Lvals]}, setup {time.time()-t0:.1f}s", flush=True)
        for a in [a_n(n), 0.75]:
            K = build_K(n, a)
            for solver in ['CLARABEL', 'SCS']:
                tv, st, Q = dec_sdp(R, G, K, solver=solver)
                print(f"  a={a:.6f} [{solver}]: max t = {tv} ({st}); 4^-n = {4.0**-n:.8e}", flush=True)
