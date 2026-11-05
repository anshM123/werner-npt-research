"""Numerical sweep (floating point; no claim rests on it): the Gershgorin quantities of G for u = 1..7, m = 1..4.

For each sample: Ghat in the frame (rho_t, c z_t, zeta'_t) (numG.Ghat), and
  A = max_{t != t'} ||Ghat_tt'||_F sigma / sqrt(s_t s_t'),     B = min_t lambda_min(Ghat_tt) / (1 - s_t/sigma)  (m >= 2),
  ratio = max_t sum_{t'} ||Ghat_tt'||_2 sqrt(s_t'/s_t) / lambda_min(Ghat_tt)   (weighted block Gershgorin, < 1 needed),
  lam = lambda_min(Ghat) (m = 1: lambda_min of the single block),
  perp = range of the eigenvalues of the zeta'-block, coup = largest entry of the zeta'-(rho, z) coupling,
  dmerge = max |Ghat(rho, z) - Ghat of the merged u = 1 system (one U-copy of weight sqrt(S))|, and the range of the
  eigenvalues of that difference (the largest one relative to the size of the merged block).
Weights: log-uniform with spreads up to e^{+-5}, plus families with S or sigma shrunk, or one U-weight dominant."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from numG import Ghat

rng = np.random.default_rng(7)


def sample(u, m, kind):
    sT = np.exp(rng.uniform(-5, 5, m))
    sU = np.exp(rng.uniform(-5, 5, u))
    if kind == 1:
        sU *= np.exp(-rng.uniform(0, 10))
    elif kind == 2:
        sT *= np.exp(-rng.uniform(0, 10))
    elif kind == 3:
        sU[0] *= np.exp(rng.uniform(0, 8))
    elif kind == 4:
        sU = np.ones(u) * np.exp(rng.uniform(-5, 5))
    return list(sT), list(sU)


t0 = time.time()
summary = {}
for u in range(1, 8):
    for m in range(1, 5):
        if u + m > 9:
            continue
        stats = dict(A=0.0, B=np.inf, ratio=0.0, lam=np.inf, perpmin=np.inf, perpmax=0.0, coup=0.0, dmerge=0.0, dmax=-np.inf, dmin=np.inf)
        N = 60 if u + m <= 7 else 25
        for it in range(N):
            sT, sU = sample(u, m, it % 5)
            sig = sum(sT)
            d = u + 1
            Gh = Ghat(sT, sU)
            stats['lam'] = min(stats['lam'], np.linalg.eigvalsh(Gh).min())
            lmin = [np.linalg.eigvalsh(Gh[t * d:(t + 1) * d, t * d:(t + 1) * d]).min() for t in range(m)]
            for t in range(m):
                row = 0.0
                for tp in range(m):
                    if tp == t:
                        continue
                    blk = Gh[t * d:(t + 1) * d, tp * d:(tp + 1) * d]
                    stats['A'] = max(stats['A'], np.linalg.norm(blk) * sig / np.sqrt(sT[t] * sT[tp]))
                    row += np.linalg.norm(blk, 2) * np.sqrt(sT[tp] / sT[t])
                if m >= 2:
                    stats['B'] = min(stats['B'], lmin[t] / (1 - sT[t] / sig))
                    stats['ratio'] = max(stats['ratio'], row / lmin[t])
            if u >= 2:
                perp = [t * d + i for t in range(m) for i in range(2, d)]
                mz = [t * d + i for t in range(m) for i in (0, 1)]
                ev = np.linalg.eigvalsh(Gh[np.ix_(perp, perp)])
                stats['perpmin'] = min(stats['perpmin'], ev.min())
                stats['perpmax'] = max(stats['perpmax'], ev.max())
                stats['coup'] = max(stats['coup'], np.abs(Gh[np.ix_(mz, perp)]).max())
                G1 = Ghat(sT, [sum(sU)])
                Dm = Gh[np.ix_(mz, mz)] - G1
                stats['dmerge'] = max(stats['dmerge'], np.abs(Dm).max())
                ev2 = np.linalg.eigvalsh((Dm + Dm.T) / 2)
                stats['dmax'] = max(stats['dmax'], ev2.max() / max(1.0, np.abs(G1).max()))
                stats['dmin'] = min(stats['dmin'], ev2.min())
        summary[(u, m)] = stats
        print(f'u={u} m={m} N={N}: A={stats["A"]:.4f} B={stats["B"]:.4f} ratio={stats["ratio"]:.4f} lam_min={stats["lam"]:.3e} '
              f'perp=[{stats["perpmin"]:.4f},{stats["perpmax"]:.4f}] coup={stats["coup"]:.2e} dmerge={stats["dmerge"]:.3f} '
              f'eig(Ghat(rho,z) - G1hat) in [{stats["dmin"]:.3f}, {stats["dmax"]:.1e}] ({time.time() - t0:.0f} s)', flush=True)
ok = all(v['lam'] > -1e-9 for v in summary.values()) and all(v['ratio'] < 1 for (u, m), v in summary.items() if m >= 2)
print('PASS' if ok else 'FAIL')
