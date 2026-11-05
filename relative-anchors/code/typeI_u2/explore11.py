"""Exploration: Gershgorin (Mc) and Schur test (Mca) with entrywise triangle bounds, exact sums (samples)."""
import numpy as np
from explore10 import entries
beta = 12 / 5
rng = np.random.default_rng(15)
wc = np.inf; wx = 0; wcomb = np.inf
for u in (2, 3, 5, 8, 12):
    for trial in range(250):
        t = rng.uniform(0, 1) if trial % 3 else 1 - np.exp(rng.uniform(-8, 0))
        if trial % 7 == 0: t = np.exp(rng.uniform(-8, 0))
        spread = rng.choice([0.3, 1, 3])
        x = np.exp(rng.uniform(-spread, spread, u))
        if trial % 4 == 1: x[0] *= np.exp(rng.uniform(1, 6))
        x = x / x.sum() * 2 * (1 - t)
        Mc, Mca, q = entries(x, t)
        w, d, rho, g, eta2, N, Lam, e, m, b, B, Pm = [q[k] for k in ('w', 'd', 'rho', 'g', 'eta2', 'N', 'Lam', 'e', 'm', 'b', 'B', 'Pm')]
        lc = np.inf
        for h in range(u):
            R = 0
            for k in range(u):
                if k == h: continue
                P1 = 0.5 - 1 / ((2 + g) * d[h] * d[k])
                P3 = 0.25 * (B + 2 * b[h] + 2 * b[k] + Pm - m[h] ** 2 / Lam[h] - m[k] ** 2 / Lam[k])
                P4 = 0.25 * (m[h] * (2 * e[h] - m[h]) / Lam[h] + m[k] * (2 * e[k] - m[k]) / Lam[k])
                R += x[k] * (P1 + P3 + abs(P4))
            lc = min(lc, Mc[h, h] - t * R)
        rows = []; cols = []
        sig = np.zeros((u, u))
        for h in range(u):
            for k in range(u):
                if h != k:
                    pos = rho[k] / (2 + g) + t * beta / (2 + g) ** 2 * eta2[k] * (eta2[h] + (N - eta2[k]) / 2)
                    neg = t / 2 * e[k] * m[k] / Lam[k]
                    sig[h, k] = max(pos, neg)
        for h in range(u):
            rows.append(abs(Mca[h, h]) + sum(x[k] * sig[h, k] for k in range(u) if k != h))
            cols.append(abs(Mca[h, h]) + sum(x[k] * sig[k, h] for k in range(u) if k != h))
        nx = np.sqrt(max(rows) * max(cols))
        wc = min(wc, lc / (1 - t)); wx = max(wx, nx / (1 - t))
        comb = lc / (1 - t) - (nx / (1 - t)) ** 2
        if comb < wcomb: wcomb = comb; arg = (t, np.round(x, 3), lc / (1 - t), nx / (1 - t))
print('Gershgorin(triangle) Mc min %.4f   Schur(triangle) Mca max %.4f   min (lc - nx^2) %.4f' % (wc, wx, wcomb), arg)
