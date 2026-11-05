"""Exploration: relaxed moment set for the T={t,t'} block. f on random convex combinations of the curve psi(y)."""
import numpy as np


def fF(al, be, F0, F1, F2, F3):
    ga = al / (2 * (1 - al)); gb = be / (2 * (1 - be))
    Gam = 1 + 0.5 * (ga + gb + F0)
    c = (1 + ga) * (1 + gb) / (2 * Gam)
    ap = 1 + al / (2 - 2 * al - be); bp = 1 + be / (2 - 2 * be - al)
    A = ap / (0.5 * (1 + ap) + 0.5 * F1); B = bp / (0.5 * (1 + bp) + 0.5 * F2)
    G3 = 1 + 0.5 * F3
    return (4 - c * c) - G3 * (A * A + B * B - c * A * B), c, A, B, G3


def psi(al, be, y):
    return np.array([1 / (2 - 2 * y), 1 / (2 - be - 2 * y), 1 / (2 - al - 2 * y), 1 / (2 - al - be - 2 * y)])


rng = np.random.default_rng(3)
worst = 1e9; which_one = 0; which_tiny = 0
for trial in range(20000):
    al, be, e = rng.dirichlet([1, 1, 1] if trial % 2 else [0.3, 0.3, 0.3])
    m = int(rng.integers(1, 6))
    ys = rng.random(m) * e
    lam = rng.dirichlet(np.ones(m))
    F = e * sum(l * psi(al, be, y) for l, y in zip(lam, ys))
    v = fF(al, be, *F)[0]
    r = v / min(1 - al, 1 - be)
    if r < worst:
        worst = r; info = (al, be, e, ys, lam)
    # which extreme is smaller
    one = fF(al, be, *(e * psi(al, be, e)))[0]; tiny = fF(al, be, *(e * psi(al, be, 0)))[0]
    if one < tiny: which_one += 1
    else: which_tiny += 1
print('relaxed hull: min f/min(1-al,1-be) =', worst, info)
print('extreme with smaller f: one copy', which_one, 'tiny copies', which_tiny)
# sign of partial derivatives 2A - cB and Gam3*A*B - 2c on random true configurations
neg1 = 0; neg2 = 0; N = 0
for trial in range(5000):
    al, be, e = rng.dirichlet([1, 1, 1])
    m = int(rng.integers(1, 6)); ys = rng.dirichlet(np.ones(m)) * e
    F = sum(y * psi(al, be, y) for y in ys)
    v, c, A, B, G3 = fF(al, be, *F)
    N += 1; neg1 += (2 * A - c * B < 0); neg2 += (G3 * A * B - 2 * c < 0)
print('2A - cB < 0 in', neg1, 'of', N, '; G3 A B - 2c < 0 in', neg2)
