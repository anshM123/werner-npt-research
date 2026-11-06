"""Exploration (floating point, no claim rests on it): growth with u of the off-diagonal constant A at general
U-weights, m = 2 (T = {t, t'} with s_t + s_t' = sigma), by random search over two-level and spread U-weights.
Reports, per u, the largest A_F = ||Ghat_tt'||_F sigma/sqrt(s_t s_t'), the spectral version A_2, and the smallest
B = lambda_min(Ghat_tt)/(1 - s_t/sigma) found, with the weights."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from numG import Ghat

rng = np.random.default_rng(7)
t0 = time.time()
for u in (3, 6, 10, 16, 24):
    best = (0, None)
    best2 = (0, None)
    bestB = (9, None)
    N = 60 if u <= 10 else 25
    for it in range(N):
        kind = it % 3
        sig = 1.0
        th = rng.uniform(0.05, 0.95)
        sT = [th * sig, (1 - th) * sig]
        S = np.exp(rng.uniform(-9, 2))
        if kind == 0:      # two values
            frac = rng.integers(1, u)
            rat = np.exp(rng.uniform(-8, 8))
            sU = np.array([rat] * frac + [1.0] * (u - frac))
        elif kind == 1:    # spread
            sU = np.exp(rng.uniform(-6, 6, u))
        else:              # one dominant, rest equal
            sU = np.ones(u)
            sU[0] = np.exp(rng.uniform(0, 8))
        sU = sU / sU.sum() * S
        d = u + 1
        Gh = Ghat(sT, list(sU))
        blk = Gh[0:d, d:2 * d]
        AF = np.linalg.norm(blk) * sig / np.sqrt(sT[0] * sT[1])
        A2 = np.linalg.norm(blk, 2) * sig / np.sqrt(sT[0] * sT[1])
        lm = np.linalg.eigvalsh(Gh[0:d, 0:d]).min() / (1 - sT[0] / sig)
        if AF > best[0]:
            best = (AF, (kind, round(th, 3), S, np.round(sU / S, 4)[:4]))
        if A2 > best2[0]:
            best2 = (A2, (kind, round(th, 3), S))
        if lm < bestB[0]:
            bestB = (lm, (kind, round(th, 3), S))
    print('u=%d: max A_F=%.4f at %s; max A_2=%.4f at %s; min B=%.4f at %s [%.0fs]' % (
        u, best[0], best[1], best2[0], best2[1], bestB[0], bestB[1], time.time() - t0), flush=True)
print('exploration done')
