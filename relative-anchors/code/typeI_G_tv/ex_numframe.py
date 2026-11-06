"""Floating-point exploration (no claim rests on it): the weighted block Gershgorin quantities of the parts of G at
two U-weight levels.  U = U1 (u1 copies of weight s1) + U2 (u2 copies of weight s2), T arbitrary.
Frame per copy t: rho_t; zeta_t = c zh_t n + zp_t nperp (+ standard parts), n = w_U/|w_U|, nperp the unit invariant
vector orthogonal to n, c = sqrt(1 + S/sigma); standard parts orthonormal.
A_inv = ||Ghat^inv_tt'||_F sigma/sqrt(s_t s_t'),  A_std_i = |Ghat^std_i_tt'| sigma/sqrt(s_t s_t'),
B_inv = lambda_min(Ghat^inv_tt)/(1 - theta_t),     B_std_i = Ghat^std_i_tt/(1 - theta_t)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from numG import Gmat


def parts(sT, u1, s1, u2, s2, scale_perp=1.0):
    sU = [s1] * u1 + [s2] * u2
    G = Gmat(sT, sU)
    m, u = len(sT), len(sU)
    d = u + 1
    wU = np.sqrt(np.array(sU))
    S = sum(sU)
    sig = sum(sT)
    S1, S2 = u1 * s1, u2 * s2
    n = wU / np.sqrt(S)
    npp = np.concatenate([np.sqrt(s1) * S2 * np.ones(u1), -np.sqrt(s2) * S1 * np.ones(u2)]) / np.sqrt(S1 * S2 * S)
    c = np.sqrt(1 + S / sig)
    P = np.zeros((d, 3))
    P[0, 0] = 1
    P[1:, 1] = c * n
    P[1:, 2] = scale_perp * npp
    f1 = np.zeros(u)
    if u1 >= 2:
        f1[0], f1[1] = 1 / np.sqrt(2), -1 / np.sqrt(2)
    f2 = np.zeros(u)
    if u2 >= 2:
        f2[u1], f2[u1 + 1] = 1 / np.sqrt(2), -1 / np.sqrt(2)
    out = dict(A_inv=[], A_std1=[], A_std2=[], B_inv=[], B_std1=[], B_std2=[])
    for t in range(m):
        for tp in range(m):
            blk = G[t * d:(t + 1) * d, tp * d:(tp + 1) * d]
            gi = P.T @ blk @ P
            g1 = f1 @ blk[1:, 1:] @ f1
            g2 = f2 @ blk[1:, 1:] @ f2
            if t == tp:
                th1 = 1 - sT[t] / sig
                out['B_inv'].append(np.linalg.eigvalsh((gi + gi.T) / 2)[0] / th1)
                if u1 >= 2:
                    out['B_std1'].append(g1 / th1)
                if u2 >= 2:
                    out['B_std2'].append(g2 / th1)
            else:
                sc = sig / np.sqrt(sT[t] * sT[tp])
                out['A_inv'].append(np.linalg.norm(gi) * sc)
                if u1 >= 2:
                    out['A_std1'].append(abs(g1) * sc)
                if u2 >= 2:
                    out['A_std2'].append(abs(g2) * sc)
    return out


if __name__ == '__main__':
    rng = np.random.default_rng(1)
    best = {}
    worst = {}
    cnt = 0
    fams = []
    for u1, u2 in [(1, 2), (1, 3), (1, 5), (2, 2), (2, 3), (3, 3), (1, 8), (4, 2), (2, 6)]:
        for rep in range(60):
            m = rng.choice([2, 2, 3])
            mode = rep % 6
            if mode == 0:
                sT = list(np.exp(rng.uniform(-4, 4, m)))
                s1, s2 = np.exp(rng.uniform(-6, 6, 2))
            elif mode == 1:     # S small vs sigma
                sT = list(np.exp(rng.uniform(-1, 1, m)))
                s1, s2 = np.exp(rng.uniform(-9, -3, 2))
            elif mode == 2:     # S large
                sT = list(np.exp(rng.uniform(-1, 1, m)))
                s1, s2 = np.exp(rng.uniform(3, 8, 2))
            elif mode == 3:     # one T-weight dominant
                sT = [1.0] + list(np.exp(rng.uniform(-7, -2, m - 1)))
                s1, s2 = np.exp(rng.uniform(-4, 4, 2))
            elif mode == 4:     # s1 tiny
                sT = list(np.exp(rng.uniform(-1, 1, m)))
                s1, s2 = np.exp(rng.uniform(-9, -4)), np.exp(rng.uniform(-2, 2))
            else:               # s1 huge
                sT = list(np.exp(rng.uniform(-1, 1, m)))
                s1, s2 = np.exp(rng.uniform(4, 9)), np.exp(rng.uniform(-2, 2))
            o = parts(sT, u1, s1, u2, s2)
            cnt += 1
            for k, v in o.items():
                if not v:
                    continue
                if k.startswith('A'):
                    if max(v) > best.get(k, (-1,))[0]:
                        best[k] = (max(v), u1, u2, [float('%.3g' % x) for x in sT], float('%.3g' % s1), float('%.3g' % s2))
                else:
                    if min(v) < worst.get(k, (1e9,))[0]:
                        worst[k] = (min(v), u1, u2, [float('%.3g' % x) for x in sT], float('%.3g' % s1), float('%.3g' % s2))
    print(cnt, 'configurations')
    for k in sorted(best):
        print('sup', k, best[k])
    for k in sorted(worst):
        print('inf', k, worst[k])
    print('exploration done')
