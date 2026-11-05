"""Exploration (floating point): the 3 x 3 block T = {t, t'} of (R_inf)_1, closed form, against j1lib, and margins."""
import numpy as np
from fractions import Fraction as Fr
from j1lib import J1, block_states


def closed_t2(x, t=0, tp=1):
    """x = s/|s| (sum 1); returns the 3x3 block in units 2^Z, order (2@t,2@t'), (2@t',2@t), (vac, 2@t 2@t')."""
    k = len(x)
    others = [u for u in range(k) if u not in (t, tp)]
    g = [xi / (2 * (1 - xi)) for xi in x]
    gam = sum(g)
    Gam = 1 + gam / 2
    a = 1 + g[t]; b = 1 + g[tp]
    ap = 1 + x[t] / (2 - 2 * x[t] - x[tp])
    bp = 1 + x[tp] / (2 - 2 * x[tp] - x[t])
    Gp = 1 + 0.5 * sum(x[u] / (2 - 2 * x[u] - x[tp]) for u in range(k) if u != tp)
    Gpp = 1 + 0.5 * sum(x[u] / (2 - 2 * x[u] - x[t]) for u in range(k) if u != t)
    G3 = 1 + 0.5 * sum(x[u] / (2 - 2 * x[u] - x[t] - x[tp]) for u in others)
    return np.array([[2, a * b / (2 * Gam), -ap / Gp], [a * b / (2 * Gam), 2, -bp / Gpp], [-ap / Gp, -bp / Gpp, 2 / G3]])


rng = np.random.default_rng(0)
for k in (4, 5, 6):
    w = rng.random(k) + 0.2
    J = J1(list(w), F=float)
    st = block_states(k, (0, 1))
    M = np.array(J.matrix(st), dtype=float)
    x = w ** 2 / np.sum(w ** 2)
    C = closed_t2(list(x))
    # scale: entries homogeneous of degree 0 in s? check ratio
    print(k, np.abs(M - C).max())
worst = 1e9
for trial in range(20000):
    k = int(rng.integers(3, 9))
    x = rng.dirichlet(np.ones(k) * rng.choice([0.1, 0.3, 1, 5]))
    x = np.maximum(x, 1e-9); x /= x.sum()
    C = closed_t2(list(x))
    ev = np.linalg.eigvalsh(C)[0]
    if ev < worst:
        worst = ev; wx = x.copy()
print('min eig', worst, np.round(wx, 4))
