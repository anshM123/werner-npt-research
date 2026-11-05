"""Exploration (floats): the virtual system (U plus one copy tau of weight sqrt(sigma)), full gauge w~^T X~ = 0.
Is the split X~ = X~' + a~ w^T (a~ perp w~, X~' w = 0, w~^T X~' = 0) F-orthogonal (as in Theorem T)?"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cf import TypeI


def null_basis(A, tol=1e-12):
    u_, s_, vt = np.linalg.svd(A)
    r = (s_ > tol * max(1, s_.max())).sum()
    return vt[r:].T


rng = np.random.default_rng(5)
for u in (2, 3, 4):
    for trial in range(3):
        wU = np.exp(rng.uniform(-1, 1, u)); wt = np.exp(rng.uniform(-2, 2))
        w = np.concatenate([[wt], wU])            # copy 0 = tau
        ti = TypeI(w, [0], F=float)
        M = np.array(ti.build())
        rows = [0] + list(range(1, u + 1))        # tau first, then U
        cols = list(range(1, u + 1))
        idx = [ti.idx[('P', r, h)] for r in rows for h in cols]
        iq = ti.idx[('Q',)]
        A = M[np.ix_(idx, idx)]; b = M[idx, iq]; Fm = A - np.outer(b, b) / M[iq, iq]
        wt_full = w[rows]; wv = w[cols]
        # coordinates: vec of (u+1) x u matrix, row-major
        def vec(X):
            return X.ravel()
        # basis of X~': X w = 0 and w~^T X = 0
        n = (u + 1) * u
        C1 = np.array([vec(np.outer(np.eye(u + 1)[r], wv)) for r in range(u + 1)])   # X w = 0 rows
        C2 = np.array([vec(np.outer(wt_full, np.eye(u)[h])) for h in range(u)])       # w~^T X = 0
        Xp = null_basis(np.vstack([C1, C2]))
        At = null_basis(wt_full[None, :])
        Om = np.array([vec(np.outer(At[:, i], wv)) for i in range(u)]).T
        cross = Xp.T @ Fm @ Om
        print(u, trial, 'dims', Xp.shape[1], Om.shape[1], ' max |cross| =', np.abs(cross).max(), ' |F| =', np.abs(Fm).max(),
              ' eig F on Omega', np.round(np.linalg.eigvalsh(Om.T @ Fm @ Om), 4), flush=True)
