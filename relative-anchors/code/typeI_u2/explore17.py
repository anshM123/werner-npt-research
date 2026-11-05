"""Exploration (floats) toward G > 0: relax the hard U-part X' of the WHOLE type I block as in Lemma S (hard part
bounded below by the Lambda-diagonal form, then minimised over X' in its subspace) and test whether the remaining
form R on (a, c, X_T, rho, zeta, q) is still >= 0 with kernel exactly the p (x) f_h directions."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cf import TypeI


def null_basis(A, tol=1e-12):
    u_, s_, vt = np.linalg.svd(A)
    r = (s_ > tol * max(1, s_.max())).sum()
    return vt[r:].T


def test(w, T):
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
    hard = [XU(np.outer(Bw[:, i], Bw[:, j])) for i in range(u - 1) for j in range(u - 1)]
    rest = [XU(np.outer(Bw[:, i], wU)) for i in range(u - 1)] + [XU(np.outer(wU, np.eye(u)[j])) for j in range(u)]
    others = [np.eye(n)[k] for k in range(n) if not (ti.labels[k][0] == 'P' and ti.labels[k][1] in U)]
    Hm = np.array(hard).T; Rm = np.array(rest + others).T
    lam = np.full((u, u), 5 / 12); np.fill_diagonal(lam, 1.5 - rho)
    # Lambda form on the hard subspace: sum lam_rh X_rh^2 with X = Bw Y Bw^T
    Lh = np.zeros((len(hard), len(hard)))
    for a_ in range(len(hard)):
        Xa = np.outer(Bw[:, a_ // (u - 1)], Bw[:, a_ % (u - 1)])
        for b_ in range(len(hard)):
            Xb = np.outer(Bw[:, b_ // (u - 1)], Bw[:, b_ % (u - 1)])
            Lh[a_, b_] = (lam * Xa * Xb).sum()
    Fh = Hm.T @ M @ Hm
    gap_h = np.linalg.eigvalsh(Fh - Lh)[0]          # Lemma H on the whole block's hard part (q = 0): >= 0
    C = Hm.T @ M @ Rm; Rr = Rm.T @ M @ Rm
    R = Rr - C.T @ np.linalg.solve(Lh, C)
    # kernel: p (x) f_h restricted to rest coordinates
    K = np.zeros((n, u))
    for j, h in enumerate(U):
        for r in range(len(w)):
            K[ti.idx[('P', r, h)], j] = w[r]
    Kr = np.linalg.lstsq(Rm, K, rcond=None)[0]
    Q, _ = np.linalg.qr(Kr)
    P = np.eye(R.shape[0]) - Q @ Q.T
    ev = np.linalg.eigvalsh(P @ R @ P + 10 * np.abs(R).max() * Q @ Q.T)
    return gap_h, ev[0], np.abs(R @ Kr).max()


if __name__ == "__main__":
    rng = np.random.default_rng(31)
    worst = np.inf; wh = np.inf
    for u in (2, 3, 4):
        for m in (1, 2, 3):
            for trial in range(25):
                spread = rng.choice([0.3, 1, 3])
                w = np.exp(rng.uniform(-spread, spread, u + m))
                if trial % 3 == 1: w[:m] *= np.exp(rng.uniform(-4, 4))
                gh, lmin, kres = test(w, list(range(m)))
                s = w ** 2
                worst = min(worst, lmin / (s[:m].sum() / s.sum())); wh = min(wh, gh)
        print(u, 'min lam(R off kernel)/(s_T/|s|) so far', round(worst, 4), ' min eig(F_hard - Lambda)', round(wh, 4), flush=True)
