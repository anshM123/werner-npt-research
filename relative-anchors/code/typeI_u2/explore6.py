"""Exploration: Phi1 (chain C1, explicit formulas, K = 1) on ALL (a, c) (no constraint a perp w) vs on a perp w.
Explicit formulas: f = sum_h[(1-t+t rho)(a^2+c^2) + 4 rho a c/(2+g)] + t/(2+g)[alpha^2+gamma^2-2p gamma - g p^2/2] - 2 p alpha/(2+g)
(valid as a formula; it equals F(Omega)/sigma only on a perp w)."""
import numpy as np


def phi1_matrix(x, t, constrained):
    x = np.asarray(x, float); u = len(x); w = np.sqrt(x); d = 1 + x; rho = x / d; g = rho.sum(); eta = w / d
    N = (eta ** 2).sum(); Lam = 1.5 - rho; e = 1 / d - eta ** 2 / (2 + g); m = rho + eta ** 2 / (2 + g)
    n = 2 * u
    def form(v):
        a, c = v[:u], v[u:]
        A, C = w * a, w * c
        p = C.sum(); al = (rho * A).sum(); ga = (rho * C).sum()
        f = ((1 - t + t * rho) * (a * a + c * c) + 4 * rho * a * c / (2 + g)).sum() \
            + t / (2 + g) * (al ** 2 + ga ** 2 - 2 * p * ga - g * p * p / 2) - 2 * p * al / (2 + g)
        G = -np.outer(eta, eta) * (2 * A[None, :] + 2 * C[:, None] + p) / (2 + g)
        dg = 2 * w * e * (a + c) - p * m
        np.fill_diagonal(G, dg)
        off = (G ** 2).sum() - (dg ** 2).sum()
        cpl = t / 4 * (12 / 5 * off + (dg ** 2 / Lam).sum())
        return f - cpl, f
    I = np.eye(n)
    M = np.zeros((n, n)); Fm = np.zeros((n, n))
    dgv = [form(I[i]) for i in range(n)]
    for i in range(n):
        M[i, i], Fm[i, i] = dgv[i]
        for j in range(i + 1, n):
            v = form(I[i] + I[j])
            M[i, j] = M[j, i] = (v[0] - dgv[i][0] - dgv[j][0]) / 2
            Fm[i, j] = Fm[j, i] = (v[1] - dgv[i][1] - dgv[j][1]) / 2
    if constrained:
        # basis of a perp w, c free
        Q, _ = np.linalg.qr(np.column_stack([w, np.eye(u)]))
        B = Q[:, 1:u]
        T = np.zeros((n, 2 * u - 1)); T[:u, :u - 1] = B; T[u:, u - 1:] = np.eye(u)
        M = T.T @ M @ T; Fm = T.T @ Fm @ T
    L = np.linalg.cholesky(Fm) if np.linalg.eigvalsh(Fm)[0] > 0 else None
    if L is None:
        return np.linalg.eigvalsh(M)[0], -1
    Li = np.linalg.inv(L)
    return np.linalg.eigvalsh(M)[0], np.linalg.eigvalsh(Li @ M @ Li.T)[0]


if __name__ == "__main__":
    rng = np.random.default_rng(2)
    for constrained in (True, False):
        worst = np.inf; worstr = np.inf
        for u in (2, 3, 5):
            for trial in range(300):
                t = rng.uniform(0, 1) if trial % 3 else 1 - np.exp(rng.uniform(-8, 0))
                spread = rng.choice([0.3, 1, 3])
                x = np.exp(rng.uniform(-spread, spread, u))
                if trial % 4 == 1: x[0] *= np.exp(rng.uniform(1, 6))
                x = x / x.sum() * 2 * (1 - t)
                lam, rel = phi1_matrix(x, t, constrained)
                worst = min(worst, lam / (1 - t)); worstr = min(worstr, rel)
        print('constrained' if constrained else 'unconstrained', 'min eig/(1-t) %.4f   min gen eig vs f %.4f' % (worst, worstr))
