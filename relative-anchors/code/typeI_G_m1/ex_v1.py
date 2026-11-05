"""Exploration (floating point): the special letter-2 direction y = w / w_t, i.e. v1 = (1, 1/s_t, ..., 1/s_t)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from ex_sep import solve
from ex_dep import Gz
np.set_printoptions(linewidth=200, precision=10, suppress=True)
rng = np.random.default_rng(11)
for u in (2, 3, 4):
    st = float(np.exp(rng.uniform(-1, 1)))
    sU = list(np.exp(rng.uniform(-1.5, 1.5, u)))
    S = sum(sU)
    v1 = np.array([1.0] + [1.0 / st] * u)
    Xs, q, X = solve(st, sU, v1)
    print(f'u={u} s_t={st:.4f} S={S:.4f} sU={np.array(sU)}')
    print(' X* scaled:\n', Xs, '\n q =', q)
    G = Gz(st, sU)
    print(' G_z v1 =', G @ v1)
    o = np.zeros(u + 1); o[1:] = 1
    print(' G_z (0,1..1) =', G @ o)
