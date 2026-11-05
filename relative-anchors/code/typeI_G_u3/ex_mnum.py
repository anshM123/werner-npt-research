"""Exploration (floating point): for m = 1, is G_z a Z-matrix, and does v = (1, beta 1) with G_z v > 0 exist?
G_z = Dg G Dg, Dg = diag(1, w_t w_h).  For each sample: max off-diagonal entry (relative), r_h = row sums of the z-block,
and the feasible interval for beta: max_h |G_{h rho}| / r_h < beta < G_{rho rho} / |sum_h G_{rho h}|."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from numG import Gmat
rng = np.random.default_rng(3)
allok = True
for u in range(1, 8):
    worst_off = -np.inf; feas = True; minr = np.inf; ratio_worst = 0.0; minv = np.inf
    for it in range(80):
        sT = [np.exp(rng.uniform(-6, 6))]
        sU = list(np.exp(rng.uniform(-6, 6, u)))
        if it % 4 == 1:
            sU = list(np.array(sU) * np.exp(-rng.uniform(0, 8)))
        if it % 4 == 2:
            sT = [sT[0] * np.exp(-rng.uniform(0, 8))]
        G = Gmat(sT, sU)
        w = np.sqrt(np.array(sT + sU))
        Dg = np.diag([1.0] + [w[0] * w[1 + h] for h in range(u)])
        Gz = Dg @ G @ Dg
        n = u + 1
        sc = np.sqrt(np.outer(np.diag(Gz), np.diag(Gz)))
        off = (Gz / sc)[~np.eye(n, dtype=bool)]
        worst_off = max(worst_off, off.max())
        r = Gz[1:, 1:].sum(axis=1)
        minr = min(minr, (r / np.diag(Gz)[1:]).min())
        lo = max(abs(Gz[1 + h, 0]) / r[h] if r[h] > 0 else np.inf for h in range(u))
        hi = Gz[0, 0] / abs(Gz[0, 1:].sum())
        feas &= lo < hi
        ratio_worst = max(ratio_worst, lo / hi)
        Ginv = np.linalg.inv(Gz)
        minv = min(minv, (Ginv / np.sqrt(np.outer(np.diag(Ginv), np.diag(Ginv)))).min())
    allok &= feas and worst_off < 0
    print(f'u = {u}: max normalised off-diagonal entry {worst_off:.3e} (<= 0: Z-matrix); min relative z-row sum {minr:.3e}; '
          f'v = (1, beta 1) feasible in all samples: {feas} (worst lo/hi = {ratio_worst:.3f}); min normalised entry of G_z^-1: {minv:.2e}', flush=True)
print('PASS' if allok else 'FAIL')
