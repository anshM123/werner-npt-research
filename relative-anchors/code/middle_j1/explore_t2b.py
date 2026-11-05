"""Exploration: det condition of the T={t,t'} block as a function of the distribution of the other weights."""
import numpy as np


def parts(al, be, ys):
    ys = np.asarray(ys, dtype=float)
    ga = al / (2 * (1 - al)); gb = be / (2 * (1 - be))
    F0 = np.sum(ys / (2 * (1 - ys))); F1 = np.sum(ys / (2 - 2 * ys - be)); F2 = np.sum(ys / (2 - 2 * ys - al))
    F3 = np.sum(ys / (2 - 2 * ys - al - be))
    Gam = 1 + 0.5 * (ga + gb + F0)
    c = (1 + ga) * (1 + gb) / (2 * Gam)
    ap = 1 + al / (2 - 2 * al - be); bp = 1 + be / (2 - 2 * be - al)
    A = ap / (0.5 * (1 + ap) + 0.5 * F1); B = bp / (0.5 * (1 + bp) + 0.5 * F2)
    G3 = 1 + 0.5 * F3
    return c, A, B, G3


def f(al, be, ys):
    c, A, B, G3 = parts(al, be, ys)
    return (4 - c * c) - G3 * (A * A + B * B - c * A * B)


rng = np.random.default_rng(1)
cnt = 0; bad = 0
for trial in range(4000):
    al, be, e = rng.dirichlet([1, 1, 1] if trial % 2 else [0.3, 0.3, 0.3])
    vals = []
    one = f(al, be, [e])
    tiny = f(al, be, [e / 2000] * 2000)
    for _ in range(5):
        m = int(rng.integers(2, 8))
        ys = rng.dirichlet(np.ones(m) * rng.choice([0.2, 1, 5])) * e
        vals.append(f(al, be, ys))
    mn = min(vals)
    cnt += 1
    if mn < min(one, tiny) - 1e-12:
        bad += 1
        if bad < 6:
            print('interior smaller', al, be, e, mn, one, tiny)
print('trials', cnt, 'interior min below both extremes:', bad)
# minimum of f over grid with one other copy and with tiny copies (normalised by (1 - be)*(1-al) ? )
worst = (1e9, None)
for al in np.linspace(0.0005, 0.999, 300):
    for be in np.linspace(0.0005, 0.999, 300):
        e = 1 - al - be
        if e < 0:
            continue
        for ys in ([e], [e / 500] * 500) if e > 0 else ([],):
            v = f(al, be, ys)
            r = v / min(1 - al, 1 - be)
            if r < worst[0]:
                worst = (r, (al, be, len(ys)))
print('min of f / min(1-al, 1-be) over grids:', worst)
