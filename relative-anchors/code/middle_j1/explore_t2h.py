"""Exploration: reduction of the T={t,t'} block to a 3-parameter check with F3 as the main variable.
Given Z = F3/e in [1/c3, 1/(c3 - 2e)], each F_i/e lies between the chord and the curve Z/(1 + d_i Z)
(d_0 = al + be, d_1 = al, d_2 = be).  f is concave in c, A, B separately, so it suffices to check the 8 endpoint
combinations."""
import numpy as np


def consts(al, be):
    ga = al / (2 * (1 - al)); gb = be / (2 * (1 - be))
    kap = (1 + ga) * (1 + gb); G0 = 1 + (ga + gb) / 2
    ap = 1 + al / (2 - 2 * al - be); bp = 1 + be / (2 - 2 * be - al)
    return kap, G0, ap, bp, (1 + ap) / 2, (1 + bp) / 2


def gmin(al, be, nl=41):
    e = 1 - al - be
    kap, G0, ap, bp, K1, K2 = consts(al, be)
    c3 = 2 - al - be
    Z0, Z1 = 1 / c3, 1 / (c3 - 2 * e)
    ds = [al + be, al, be]
    out = 1e9
    for lam in np.linspace(0, 1, nl):
        Z = Z0 + lam * (Z1 - Z0)
        rng_ = []
        for d in ds:
            cur = Z / (1 + d * Z)
            x0, x1 = Z0 / (1 + d * Z0), Z1 / (1 + d * Z1)
            ch = x0 + (x1 - x0) * lam
            rng_.append((e * ch, e * cur))
        G3 = 1 + e * Z / 2
        for F0 in rng_[0]:
            c = kap / (2 * (G0 + F0 / 2))
            for F1 in rng_[1]:
                A = ap / (K1 + F1 / 2)
                for F2 in rng_[2]:
                    B = bp / (K2 + F2 / 2)
                    out = min(out, 4 - c * c - G3 * (A * A + B * B - c * A * B))
    return out


rng = np.random.default_rng(11)
worst = 1e9
for al in np.linspace(1e-4, 1 - 2e-4, 150):
    for be in np.linspace(1e-4, 1 - 2e-4, 150):
        if al + be >= 1:
            continue
        v = gmin(al, be) / min(1 - al, 1 - be)
        if v < worst:
            worst = v; info = (al, be)
print('grid: min g / min(1-al,1-be) =', worst, info)
for trial in range(5000):
    al, be, e = rng.dirichlet([0.2, 0.2, 0.2])
    if e < 1e-12: continue
    v = gmin(al, be) / min(1 - al, 1 - be)
    if v < worst:
        worst = v; info = (al, be)
print('random: min g / min(1-al,1-be) =', worst, info)
