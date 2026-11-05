"""Exploration: the reduced 3-parameter check for the T={t,t'} block.
g(c, B, F1) = 4 - c^2 - (1 + chord31(F1)/2)(A(F1)^2 + B^2 - c A(F1) B), c in {c_lo, c_hi}, B in {B_lo, B_hi}."""
import numpy as np


def consts(al, be):
    ga = al / (2 * (1 - al)); gb = be / (2 * (1 - be))
    kap = (1 + ga) * (1 + gb); G0 = 1 + (ga + gb) / 2
    ap = 1 + al / (2 - 2 * al - be); bp = 1 + be / (2 - 2 * be - al)
    return kap, G0, ap, bp, (1 + ap) / 2, (1 + bp) / 2


def gmin(al, be, nl=41):
    e = 1 - al - be
    kap, G0, ap, bp, K1, K2 = consts(al, be)
    c0, c1, c2, c3 = 2, 2 - be, 2 - al, 2 - al - be
    F0r = (e / c0, e / (c0 - 2 * e)); F2r = (e / c2, e / (c2 - 2 * e)); F1r = (e / c1, e / (c1 - 2 * e))
    F3r = (e / c3, e / (c3 - 2 * e))
    cs = [kap / (2 * (G0 + F0 / 2)) for F0 in F0r]
    Bs = [bp / (K2 + F2 / 2) for F2 in F2r]
    out = 1e9
    for lam in np.linspace(0, 1, nl):
        F1 = F1r[0] + lam * (F1r[1] - F1r[0])
        F3 = F3r[0] + lam * (F3r[1] - F3r[0])  # chord in the (F1, F3) plane
        A = ap / (K1 + F1 / 2); G3 = 1 + F3 / 2
        for c in cs:
            for B in Bs:
                out = min(out, 4 - c * c - G3 * (A * A + B * B - c * A * B))
    return out


rng = np.random.default_rng(11)
worst = 1e9
for al in np.linspace(1e-4, 1 - 2e-4, 120):
    for be in np.linspace(1e-4, 1 - 2e-4, 120):
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
