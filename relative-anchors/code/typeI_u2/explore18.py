"""Exploration (floats): as explore17, but the hard part is split X' = X'' + tau P/(u-1) (P = I - w w^T/S, trace
direction kept exactly as a 'rest' variable) and only the trace-free X'' is relaxed (Lambda form on X'')."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cf import TypeI
from explore17 import null_basis


def test(w, T, keep_diag=False):
    ti = TypeI(w, T, F=float)
    M = np.array(ti.build()); n = ti.n
    U = ti.U; u = len(U); w = np.asarray(w, float); s = w ** 2
    wU = w[U]; S = s[U].sum(); sig = s[T].sum()
    d = s[U] + S / 2 + sig; rho = s[U] / d
    Bw = null_basis(wU[None, :])
    def XU(X):
        v = np.zeros(n)
        for i, r in enumerate(U):
            for j, h in enumerate(U):
                v[ti.idx[('P', r, h)]] = X[i, j]
        return v
    Pm = np.eye(u) - np.outer(wU, wU) / S
    # basis of V' (dim (u-1)^2) orthogonal to P: trace-free part
    Vb = [np.outer(Bw[:, i], Bw[:, j]) for i in range(u - 1) for j in range(u - 1)]
    Vmat = np.array([V.ravel() for V in Vb]).T
    pn = Pm.ravel() / np.linalg.norm(Pm)
    Vperp = Vmat - np.outer(pn, pn @ Vmat)
    Qv, Rv = np.linalg.qr(Vperp)
    rank = (np.abs(np.diag(Rv)) > 1e-10).sum()
    Hbasis = [Qv[:, i].reshape(u, u) for i in range(rank)]
    hard = [XU(H) for H in Hbasis]
    rest = [XU(Pm / (u - 1))] + [XU(np.outer(Bw[:, i], wU)) for i in range(u - 1)] + [XU(np.outer(wU, np.eye(u)[j])) for j in range(u)]
    others = [np.eye(n)[k] for k in range(n) if not (ti.labels[k][0] == 'P' and ti.labels[k][1] in U)]
    Hm = np.array(hard).T; Rm = np.array(rest + others).T
    lam = np.full((u, u), 5 / 12); np.fill_diagonal(lam, 1.5 - rho)
    Lh = np.array([[(lam * A * B).sum() for B in Hbasis] for A in Hbasis])
    Rr = Rm.T @ M @ Rm
    if rank == 0:
        R = Rr
    else:
        C = Hm.T @ M @ Rm
        R = Rr - C.T @ np.linalg.solve(Lh, C)
    K = np.zeros((n, u))
    for j, h in enumerate(U):
        for r in range(len(w)):
            K[ti.idx[('P', r, h)], j] = w[r]
    Kr = np.linalg.lstsq(Rm, K, rcond=None)[0]
    Q, _ = np.linalg.qr(Kr)
    P = np.eye(R.shape[0]) - Q @ Q.T
    ev = np.linalg.eigvalsh(P @ R @ P + 10 * np.abs(R).max() * Q @ Q.T)
    return ev[0]


rng = np.random.default_rng(31)
for u in (2, 3, 4, 5):
    worst = np.inf
    for m in (1, 2, 3):
        for trial in range(25):
            spread = rng.choice([0.3, 1, 3])
            w = np.exp(rng.uniform(-spread, spread, u + m))
            if trial % 3 == 1: w[:m] *= np.exp(rng.uniform(-4, 4))
            lmin = test(w, list(range(m)))
            s = w ** 2
            worst = min(worst, lmin / (s[:m].sum() / s.sum()))
    print(u, 'trace kept: min lam(R off kernel)/(s_T/|s|)', round(worst, 4), flush=True)
