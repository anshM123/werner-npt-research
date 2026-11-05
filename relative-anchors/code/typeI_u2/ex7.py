"""Where is inf lambda_min(C_u) / (s_T/|s|) = 1/2 approached?  Targeted search."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from sweep import analyse
rng = np.random.default_rng(3)
best = (9, None)
for trial in range(3000):
    u = int(rng.integers(2, 5)); m = int(rng.integers(1, 3)); k = u + m
    w = np.exp(rng.uniform(-4, 4, k) * rng.choice([0.2, 1, 2]))
    r = analyse(w, list(range(m)))
    if r['relC'] < best[0]:
        s = w ** 2
        best = (r['relC'], (u, m, np.round(s / s.sum(), 5)))
print("min ratio", best)
# scan: s_T -> 0 or s_T -> large, one U weight dominating
for sT in (1e-6, 1e-3, 1e-1, 1, 10, 1e3, 1e6):
    for wU in ([1, 1], [1, 1e-3], [1, 1, 1], [1, 1e-3, 1e-3]):
        w = np.array([np.sqrt(sT)] + list(wU), float)
        r = analyse(w, [0])
        print(f"s_T={sT:g} w_U={wU}: lam(C)/(s_T/|s|) = {r['relC']:.4f}")
