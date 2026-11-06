"""Floating-point sweep (no claim rests on it): the weighted block Gershgorin quantities of Proposition U1 / Theorem U2
for GENERAL U-weights, u = 3..10, m = 2, 3 (k up to 13), in the frame of numG.Ghat (rho_t, c z_t, zeta'_t), c = sqrt(1 + S/sigma):
  A = max_{t != t'} ||Ghat_tt'||_F sigma/sqrt(s_t s_t'),   B = min_t lambda_min(Ghat_tt)/(1 - s_t/sigma),
  ratio = max_t sum_{t'} ||Ghat_tt'||_2 sqrt(s_t'/s_t)/lambda_min(Ghat_tt).
Families: log-uniform weights (spread e^{+-5}); S shrunk; sigma shrunk; one U-weight dominant; one U-weight tiny; two
U-weight values; equal U-weights.  Reports the worst values and where they occur.  The run passes if every sample
has sup A < 1/2 < inf B (the window of Theorem EM) and ratio < 1."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from numG import Ghat

rng = np.random.default_rng(20261105)


def sample(u, m, kind):
    sT = np.exp(rng.uniform(-5, 5, m))
    sU = np.exp(rng.uniform(-5, 5, u))
    if kind == 1:
        sU *= np.exp(-rng.uniform(4, 10))
    elif kind == 2:
        sT *= np.exp(-rng.uniform(4, 10))
    elif kind == 3:
        sU[0] *= np.exp(rng.uniform(4, 8))
    elif kind == 4:
        sU[0] *= np.exp(-rng.uniform(4, 10))
    elif kind == 5:
        v = np.exp(rng.uniform(-5, 5, 2))
        sU = np.where(np.arange(u) < u // 2, v[0], v[1])
    elif kind == 6:
        sU = np.ones(u) * np.exp(rng.uniform(-5, 5))
        sU *= np.exp(-rng.uniform(0, 8))
    return list(sT), list(sU)


t0 = time.time()
worst = dict(A=0.0, B=np.inf, ratio=0.0, lam=np.inf)
where = {}
ok = True
for u in (3, 4, 5, 6, 8, 10):
    for m in (2, 3):
        N = 28 if u <= 6 else 14
        for it in range(N):
            kind = it % 7
            sT, sU = sample(u, m, kind)
            sig = sum(sT)
            d = u + 1
            Gh = Ghat(sT, sU)
            lam = np.linalg.eigvalsh(Gh).min()
            lmin = [np.linalg.eigvalsh(Gh[t * d:(t + 1) * d, t * d:(t + 1) * d]).min() for t in range(m)]
            for t in range(m):
                row = 0.0
                for tp in range(m):
                    if tp == t:
                        continue
                    blk = Gh[t * d:(t + 1) * d, tp * d:(tp + 1) * d]
                    A = np.linalg.norm(blk) * sig / np.sqrt(sT[t] * sT[tp])
                    if A > worst['A']:
                        worst['A'] = A
                        where['A'] = (u, m, kind)
                    row += np.linalg.norm(blk, 2) * np.sqrt(sT[tp] / sT[t])
                B = lmin[t] / (1 - sT[t] / sig)
                if B < worst['B']:
                    worst['B'] = B
                    where['B'] = (u, m, kind)
                worst['ratio'] = max(worst['ratio'], row / lmin[t])
            worst['lam'] = min(worst['lam'], lam)
        print('u=%d m=%d: worst so far A=%.4f (at %s) B=%.4f (at %s) ratio=%.4f lam_min=%.2e [%.0fs]' % (
            u, m, worst['A'], where.get('A'), worst['B'], where.get('B'), worst['ratio'], worst['lam'],
            time.time() - t0), flush=True)
ok = worst['A'] < 0.5 < worst['B'] and worst['ratio'] < 1 and worst['lam'] > -1e-9
print('families: 0 log-uniform, 1 S shrunk, 2 sigma shrunk, 3 one U-weight dominant, 4 one U-weight tiny, 5 two values, '
      '6 equal U-weights')
print('PASS' if ok else 'FAIL')
