"""Exploration: margin of chain C1 as a function of t = sigma/K (K = S/2 + sigma)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from explore3 import mats
rng = np.random.default_rng(9)
for t in (0.02, 0.1, 0.3, 0.5, 0.7, 0.9, 0.97, 0.995):
    worst = np.inf
    for u in (2, 3, 4, 6):
        for trial in range(120):
            spread = rng.choice([0.3, 1, 2, 4])
            w = np.exp(rng.uniform(-spread, spread, u))
            if trial % 4 == 1:
                w[0] *= np.exp(rng.uniform(2, 5))
            S = (w ** 2).sum()
            # t = sig/(S/2 + sig) -> sig = t S / (2 (1 - t))
            sig = t * S / (2 * (1 - t))
            worst = min(worst, mats(w, sig)[0])
    print('t = %.3f  min Phi1/f = %.4f' % (t, worst), flush=True)
