import numpy as np
from explore7 import blocks
rng = np.random.default_rng(12)
mins = dict(la=np.inf, lc=np.inf); mx = 0
for u in (2, 3, 5, 8):
    for trial in range(300):
        t = rng.uniform(0, 1) if trial % 3 else 1 - np.exp(rng.uniform(-8, 0))
        if trial % 7 == 0: t = np.exp(rng.uniform(-8, 0))
        spread = rng.choice([0.3, 1, 3])
        x = np.exp(rng.uniform(-spread, spread, u))
        if trial % 4 == 1: x[0] *= np.exp(rng.uniform(1, 6))
        x = x / x.sum() * 2 * (1 - t)
        Ma, Mc, Mca = blocks(x, t)
        la, lc = np.linalg.eigvalsh(Ma)[0] / (1 - t), np.linalg.eigvalsh(Mc)[0] / (1 - t)
        nc = np.linalg.norm(Mca, 2) / (1 - t)
        if la < mins['la']: mins['la'] = la; aa = (t, np.round(x, 3))
        if lc < mins['lc']: mins['lc'] = lc; ac = (t, np.round(x, 3))
        if nc > mx: mx = nc; ax = (t, np.round(x, 3))
print('min lam(Ma)/(1-t)', mins['la'], aa)
print('min lam(Mc)/(1-t)', mins['lc'], ac)
print('max ||Mca||/(1-t)', mx, ax)
