"""Exploration: weighted Gershgorin / Schur tests on the exact Phi1 blocks (weights w_h).
 Mc: min_h [Mc_hh - sum_{k!=h} |Mc_hk| w_k / w_h] / (1-t)
 Mca: weighted Schur test bound sqrt(max row sum * max col sum) / (1-t) with weights w."""
import numpy as np
from explore7 import blocks
rng = np.random.default_rng(13)
worst_c = np.inf; worst_x = 0; worst_comb = np.inf
for u in (2, 3, 5, 8, 12):
    for trial in range(250):
        t = rng.uniform(0, 1) if trial % 3 else 1 - np.exp(rng.uniform(-8, 0))
        if trial % 7 == 0: t = np.exp(rng.uniform(-8, 0))
        spread = rng.choice([0.3, 1, 3])
        x = np.exp(rng.uniform(-spread, spread, u))
        if trial % 4 == 1: x[0] *= np.exp(rng.uniform(1, 6))
        x = x / x.sum() * 2 * (1 - t)
        w = np.sqrt(x)
        Ma, Mc, Mca = blocks(x, t)
        off = np.abs(Mc - np.diag(np.diag(Mc)))
        gc = min(Mc[h, h] - (off[h] * w).sum() / w[h] for h in range(u)) / (1 - t)
        A = np.abs(Mca)
        row = max((A[h] * w).sum() / w[h] for h in range(u))
        col = max((A[:, k] * w).sum() / w[k] for k in range(u))
        nx = np.sqrt(row * col) / (1 - t)
        # unweighted Schur test too
        row1 = A.sum(1).max(); col1 = A.sum(0).max(); nx1 = np.sqrt(row1 * col1) / (1 - t)
        nx = min(nx, nx1)
        worst_c = min(worst_c, gc); worst_x = max(worst_x, nx)
        la = np.linalg.eigvalsh(Ma)[0] / (1 - t)
        worst_comb = min(worst_comb, la * gc - nx ** 2)
        if la * gc - nx ** 2 <= worst_comb: arg = (t, np.round(x, 3), la, gc, nx)
print('Gershgorin Mc min:', worst_c, '  Schur-test Mca max:', worst_x, '  min la*gc - nx^2:', worst_comb, arg)
