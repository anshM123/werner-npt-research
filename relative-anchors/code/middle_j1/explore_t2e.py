"""Exploration: is Phi(F0,F1,F2) = (4 - c^2)/(A^2 + B^2 - cAB) - 1 concave on the relevant region?"""
import numpy as np


def consts(al, be):
    ga = al / (2 * (1 - al)); gb = be / (2 * (1 - be))
    kap = (1 + ga) * (1 + gb); G0 = 1 + (ga + gb) / 2
    ap = 1 + al / (2 - 2 * al - be); bp = 1 + be / (2 - 2 * be - al)
    return kap, G0, ap, bp, (1 + ap) / 2, (1 + bp) / 2


def Phi(al, be, F):
    kap, G0, ap, bp, K1, K2 = consts(al, be)
    F0, F1, F2 = F
    c = kap / (2 * (G0 + F0 / 2)); A = ap / (K1 + F1 / 2); B = bp / (K2 + F2 / 2)
    return (4 - c * c) / (A * A + B * B - c * A * B) - 1


rng = np.random.default_rng(5)
worst = -1e9; nbad = 0; N = 0
for trial in range(20000):
    al, be, e = rng.dirichlet([1, 1, 1] if trial % 2 else [0.3, 0.3, 0.3])
    # random point in the relevant region: from a random configuration
    m = int(rng.integers(1, 6)); ys = rng.dirichlet(np.ones(m)) * e
    F = np.array([np.sum(ys / (2 - 2 * ys)), np.sum(ys / (2 - be - 2 * ys)), np.sum(ys / (2 - al - 2 * ys))])
    h = 1e-4 * max(1e-3, F.max())
    H = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            ei = np.eye(3)[i] * h; ej = np.eye(3)[j] * h
            H[i, j] = (Phi(al, be, F + ei + ej) - Phi(al, be, F + ei - ej) - Phi(al, be, F - ei + ej) + Phi(al, be, F - ei - ej)) / (4 * h * h)
    ev = np.linalg.eigvalsh(H)[-1]
    N += 1
    scale = np.abs(H).max() + 1e-12
    if ev > 1e-6 * scale:
        nbad += 1
        if ev / scale > worst:
            worst = ev / scale; info = (al, be, e, F, np.linalg.eigvalsh(H))
print('points', N, 'with a positive Hessian eigenvalue:', nbad)
if nbad:
    print('worst relative', worst, info)
