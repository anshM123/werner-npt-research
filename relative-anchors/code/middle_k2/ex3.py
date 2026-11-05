"""Exploration: Schur complement on the soft space S (R eliminated), rescaled by 1/sqrt(x_h).
Also: which states carry the R-component of the soft eigenvectors (first order in w_h)."""
import numpy as np, sys
from ex1 import reduced, kernel_basis
from ex2 import soft_basis
np.set_printoptions(linewidth=220, precision=4, suppress=True)


def schur_soft(w):
    m, Y, F = reduced(list(w))
    k, n = m.k, m.np
    s = np.array(w, float) ** 2; x = s / s.sum()
    K = kernel_basis(m)
    Sb = [np.linalg.qr(B)[0] for B in soft_basis(m)]
    allS = np.hstack(Sb)
    Q, _ = np.linalg.qr(np.hstack([K, allS]))
    Rb = np.linalg.svd(np.eye(k * n) - Q @ Q.T)[0][:, :k * n - Q.shape[1]]
    FRR = Rb.T @ F @ Rb; FSR = allS.T @ F @ Rb; FSS = allS.T @ F @ allS
    Sig = FSS - FSR @ np.linalg.solve(FRR, FSR.T)
    d = np.repeat(1 / np.sqrt(x), k - 2)
    Sh = d[:, None] * Sig * d[None, :]
    return m, x, Sh, np.linalg.eigvalsh(FRR)[0]


if __name__ == '__main__':
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    worst = {}
    for k in (4, 5, 6):
        res = []
        for trial in range(60 if k < 6 else 15):
            mode = trial % 5
            if mode == 0:
                w = rng.random(k) + 0.3
            elif mode == 1:
                w = np.exp(1.5 * rng.normal(size=k))
            elif mode == 2:
                w = np.exp(rng.normal(size=k)); w[0] *= 1e-2
            elif mode == 3:
                w = np.exp(rng.normal(size=k)); w[:2] *= 1e-2
            else:
                w = np.exp(rng.normal(size=k)); w[0] *= 30
            m, x, Sh, lamR = schur_soft(w)
            b = k - 2
            lam = np.linalg.eigvalsh(Sh)[0]
            dmin = min(np.linalg.eigvalsh(Sh[h * b:(h + 1) * b, h * b:(h + 1) * b])[0] for h in range(k))
            # weighted Gershgorin ratio with weights sqrt(x)
            ratio = 0
            for h in range(k):
                lh = np.linalg.eigvalsh(Sh[h * b:(h + 1) * b, h * b:(h + 1) * b])[0]
                off = sum(np.linalg.norm(Sh[h * b:(h + 1) * b, g * b:(g + 1) * b], 2) * np.sqrt(x[g] / x[h]) for g in range(k) if g != h)
                ratio = max(ratio, off / lh)
            # plain (unweighted) block Gershgorin
            ratio2 = 0
            for h in range(k):
                lh = np.linalg.eigvalsh(Sh[h * b:(h + 1) * b, h * b:(h + 1) * b])[0]
                off = sum(np.linalg.norm(Sh[h * b:(h + 1) * b, g * b:(g + 1) * b], 2) for g in range(k) if g != h)
                ratio2 = max(ratio2, off / lh)
            offn = max(np.linalg.norm(Sh[h * b:(h + 1) * b, g * b:(g + 1) * b], 2) for h in range(k) for g in range(k) if g != h)
            res.append((lam, dmin, ratio, ratio2, lamR, offn, x.min(), x.max()))
        res = np.array(res)
        print(f"k={k}: min lambda(Sigma_hat) {res[:, 0].min():.4f}; min diag-block eig {res[:, 1].min():.4f}; "
              f"max weighted-Gersh ratio {res[:, 2].max():.4f}; max plain ratio {res[:, 3].max():.4f}; min lambda(F_RR) {res[:, 4].min():.4f}; "
              f"max off-block norm {res[:, 5].max():.4f}")
