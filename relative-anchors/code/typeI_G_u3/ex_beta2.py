"""Exploration (floating point): m = 1, the window [lo, hi] for beta, in units of 1/s_t, against x = S/s_t."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from numG import Gmat
for u, shape in ((1, None), (2, 'eq'), (3, 'eq'), (3, 'spread'), (5, 'spread')):
    print(f'u = {u}, U-weights {shape}')
    for x in (1e-4, 1e-2, 0.1, 0.5, 1, 2, 10, 100, 1e4):
        st = 1.0
        if shape == 'spread':
            base = np.exp(np.linspace(-3, 3, u))
        else:
            base = np.ones(u)
        sU = x * base / base.sum()
        G = Gmat([st], list(sU))
        w = np.sqrt(np.array([st] + list(sU)))
        Dg = np.diag([1.0] + [w[0] * w[1 + h] for h in range(u)])
        Gz = Dg @ G @ Dg
        r = Gz[1:, 1:].sum(axis=1)
        lo = max(abs(Gz[1 + h, 0]) / r[h] for h in range(u))
        hi = Gz[0, 0] / abs(Gz[0, 1:].sum())
        print(f'   S/s_t = {x:8.0e}: lo*s_t = {lo:.4f}, hi*s_t = {hi:.4f}, sqrt(lo hi)*s_t = {np.sqrt(lo*hi):.4f}')
