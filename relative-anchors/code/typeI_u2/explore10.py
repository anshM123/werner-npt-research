"""Exploration: explicit entry formulas for Mc, Mca (K = 1) checked against blocks(); then Gershgorin / Schur
bounds with entrywise triangle bounds |pi_hk| <= max(P1+P3, P4) etc."""
import numpy as np
from explore7 import blocks
beta = 12 / 5


def entries(x, t):
    x = np.asarray(x, float); u = len(x); w = np.sqrt(x); d = 1 + x; rho = x / d; g = rho.sum(); eta2 = x / d ** 2
    N = eta2.sum(); Lam = 1.5 - rho; e = 1 / d - eta2 / (2 + g); m = rho + eta2 / (2 + g)
    b = beta * eta2 * (N - eta2) / (2 + g) ** 2; B = b.sum(); Pm = (m ** 2 / Lam).sum()
    Mc = np.zeros((u, u)); Mca = np.zeros((u, u))
    for h in range(u):
        Mc[h, h] = (1 - t) + t * rho[h] + t * eta2[h] / (2 + g) - t * x[h] / 2 - t / 4 * x[h] * (8 * b[h] + B + 4 * e[h] * (e[h] - m[h]) / Lam[h] + Pm)
        for k in range(u):
            if k != h:
                Mc[h, k] = w[h] * w[k] * t * (1 / ((2 + g) * d[h] * d[k]) - 0.5 - 0.25 * (B + 2 * b[h] + 2 * b[k] + Pm - 2 * e[h] * m[h] / Lam[h] - 2 * e[k] * m[k] / Lam[k]))
    kap = w * (e * m / Lam - beta * eta2 * (N - eta2) / (2 + g) ** 2)
    kh = t * x * e ** 2 / Lam - t * beta / (2 + g) ** 2 * eta2 ** 2 * x
    zeta = eta2 * w
    for h in range(u):
        for k in range(u):
            Mca[h, k] = (2 * (h == k) - w[h] * w[k]) * rho[k] / (2 + g) - (h == k) * kh[h] - t * beta / (2 + g) ** 2 * zeta[h] * zeta[k] + t / 2 * w[h] * kap[k]
    return Mc, Mca, dict(w=w, d=d, rho=rho, g=g, eta2=eta2, N=N, Lam=Lam, e=e, m=m, b=b, B=B, Pm=Pm)


if __name__ == "__main__":
  rng = np.random.default_rng(14)
  err = 0; worst = np.inf
  for u in (2, 3, 5, 8):
      for trial in range(120):
          t = rng.uniform(0, 1) if trial % 3 else 1 - np.exp(rng.uniform(-8, 0))
          spread = rng.choice([0.3, 1, 3])
          x = np.exp(rng.uniform(-spread, spread, u))
          if trial % 4 == 1: x[0] *= np.exp(rng.uniform(1, 6))
          x = x / x.sum() * 2 * (1 - t)
          Ma, Mc0, Mca0 = blocks(x, t)
          Mc, Mca, q = entries(x, t)
          err = max(err, np.abs(Mc - Mc0).max() / (1 - t), np.abs(Mca - Mca0).max() / (1 - t))
  print('max relative entry error', err)
