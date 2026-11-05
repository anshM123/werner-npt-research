"""Exploration: relaxations of the moment constraints for the T={t,t'} block.
R_B: F0, F1, F2 free in their ranges [e/c_i, e/(c_i - 2e)], F3 = min(chord31(F1), chord32(F2))."""
import numpy as np


def consts(al, be):
    ga = al / (2 * (1 - al)); gb = be / (2 * (1 - be))
    kap = (1 + ga) * (1 + gb); G0 = 1 + (ga + gb) / 2
    ap = 1 + al / (2 - 2 * al - be); bp = 1 + be / (2 - 2 * be - al)
    return kap, G0, ap, bp, (1 + ap) / 2, (1 + bp) / 2


def f(al, be, F0, F1, F2, F3):
    kap, G0, ap, bp, K1, K2 = consts(al, be)
    c = kap / (2 * (G0 + F0 / 2)); A = ap / (K1 + F1 / 2); B = bp / (K2 + F2 / 2); G3 = 1 + F3 / 2
    return 4 - c * c - G3 * (A * A + B * B - c * A * B)


def chord(e, ci, cj, Fi):
    """upper chord bound of F_j in terms of F_i (j more singular): through tiny (e/ci, e/cj) and one copy."""
    x0, z0 = e / ci, e / cj
    x1, z1 = e / (ci - 2 * e), e / (cj - 2 * e)
    return z0 + (z1 - z0) / (x1 - x0) * (Fi - x0)


rng = np.random.default_rng(7)
worst = 1e9
for trial in range(4000):
    al, be, e = rng.dirichlet([1, 1, 1] if trial % 2 else [0.3, 0.3, 0.3])
    c0, c1, c2, c3 = 2, 2 - be, 2 - al, 2 - al - be
    for F0 in np.linspace(e / c0, e / (c0 - 2 * e), 5):
        for F1 in np.linspace(e / c1, e / (c1 - 2 * e), 7):
            for F2 in np.linspace(e / c2, e / (c2 - 2 * e), 7):
                F3 = min(chord(e, c1, c3, F1), chord(e, c2, c3, F2))
                v = f(al, be, F0, F1, F2, F3) / min(1 - al, 1 - be)
                if v < worst:
                    worst = v; info = (al, be, e, F0, F1, F2, F3)
print('R_B: min f/min(1-al,1-be) =', worst, info)
