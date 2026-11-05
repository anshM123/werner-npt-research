import numpy as np
from rowbounds2 import bounds
from explore10 import entries
rng = np.random.default_rng(16)
for u in (2, 3, 5, 8):
    for trial in range(200):
        t = rng.uniform(0, 1) if trial % 3 else 1 - np.exp(rng.uniform(-8, 0))
        spread = rng.choice([0.3, 1, 3])
        x = np.exp(rng.uniform(-spread, spread, u))
        if trial % 4 == 1: x[0] *= np.exp(rng.uniform(1, 6))
        S = 2 * (1 - t); x = x / x.sum() * S
        Mc, Mca, q = entries(x, t)
        w = q['w']
        for h in range(u):
            ex = (Mc[h, h] + sum(Mc[h, k] * w[k] / w[h] for k in range(u) if k != h)) / (1 - t)
            Sp = S - x[h]; rSp = Sp / (1 + Sp); gp = q['g'] - x[h] / (1 + x[h])
            theta = (gp - rSp) / (Sp - rSp) if Sp - rSp > 1e-15 else 0.0
            Nub_p = gp - (Sp - gp) / (1 + Sp)
            nu = (q['N'] - q['eta2'][h]) / Nub_p if Nub_p > 0 else 0.0
            b = bounds(S, x[h] / S, theta, nu, one=1.0)
            rowex = sum(abs(Mca[h, k]) * w[k] / w[h] for k in range(u)) / (1 - t)
            colex = sum(abs(Mca[k, h]) * w[k] / w[h] for k in range(u)) / (1 - t)
            bad = []
            if b['LB'] > ex + 1e-9: bad.append(('LB', b['LB'], ex))
            if max(b['ROWp'], b['ROWm']) < rowex - 1e-9: bad.append(('ROW', max(b['ROWp'], b['ROWm']), rowex))
            if max(b['COLp'], b['COLm']) < colex - 1e-9: bad.append(('COL', max(b['COLp'], b['COLm']), colex))
            if not (-1e-9 <= nu <= 1 + 1e-9 and -1e-9 <= theta <= 1 + 1e-9): bad.append(('range', theta, nu))
            if bad: print(u, t, np.round(x, 5), h, bad)
