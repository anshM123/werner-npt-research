"""Exploration: blocks of Phi1 (unconstrained, K = 1): Phi1 = a^T Ma a + c^T Mc c + 2 c^T Mca a.
Check the 'three norms' sufficient condition  lam_min(Ma) lam_min(Mc) >= ||Mca||^2, and record worst ratios."""
import numpy as np
from explore6 import phi1_matrix


def blocks(x, t):
    x = np.asarray(x, float); u = len(x); w = np.sqrt(x); d = 1 + x; rho = x / d; g = rho.sum(); eta = w / d
    Lam = 1.5 - rho; e = 1 / d - eta ** 2 / (2 + g); m = rho + eta ** 2 / (2 + g)
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
        return f - t / 4 * (12 / 5 * off + (dg ** 2 / Lam).sum())
    I = np.eye(n); M = np.zeros((n, n))
    dv = [form(I[i]) for i in range(n)]
    for i in range(n):
        M[i, i] = dv[i]
        for j in range(i + 1, n):
            M[i, j] = M[j, i] = (form(I[i] + I[j]) - dv[i] - dv[j]) / 2
    return M[:u, :u], M[u:, u:], M[u:, :u]


if __name__ == "__main__":
  rng = np.random.default_rng(11)
  worst = np.inf
  for u in (2, 3, 5, 8):
      for trial in range(300):
          t = rng.uniform(0, 1) if trial % 3 else 1 - np.exp(rng.uniform(-8, 0))
          spread = rng.choice([0.3, 1, 3])
          x = np.exp(rng.uniform(-spread, spread, u))
          if trial % 4 == 1: x[0] *= np.exp(rng.uniform(1, 6))
          x = x / x.sum() * 2 * (1 - t)
          Ma, Mc, Mca = blocks(x, t)
          la, lc = np.linalg.eigvalsh(Ma)[0], np.linalg.eigvalsh(Mc)[0]
          nc = np.linalg.norm(Mca, 2)
          r = np.sqrt(max(la, 0) * max(lc, 0)) / nc if nc > 0 else np.inf
          if r < worst:
              worst = r; arg = (t, np.round(x, 3), la / (1 - t), lc / (1 - t), nc / (1 - t))
  print('worst sqrt(la lc)/||Mca|| =', worst, 'at', arg)
