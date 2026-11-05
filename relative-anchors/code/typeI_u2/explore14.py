import numpy as np
from rowbounds2 import bounds, pair_Q
from explore10 import entries
# (1) validity on samples: bound <= exact signed Gershgorin row; pi <= 0; ROW/COL bounds >= exact
rng = np.random.default_rng(16)
viol = 0; maxpi = -np.inf; worst_ex = np.inf
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
            pis = [Mc[h, k] / (w[h] * w[k] * t) for k in range(u) if k != h]
            maxpi = max(maxpi, max(pis))
            ex = (Mc[h, h] + sum(Mc[h, k] * w[k] / w[h] for k in range(u) if k != h)) / (1 - t)
            worst_ex = min(worst_ex, ex)
            Sp = S - x[h]; rSp = Sp / (1 + Sp); gp = q['g'] - x[h] / (1 + x[h])
            theta = (gp - rSp) / (Sp - rSp) if Sp - rSp > 1e-15 else 0.0
            Nub_p = gp - (Sp - gp) / (1 + Sp)
            nu = (q['N'] - q['eta2'][h]) / Nub_p if Nub_p > 0 else 0.0
            b = bounds(S, x[h] / S, theta, nu, one=1.0)
            if b['LB'] > ex + 1e-9: viol += 1
            rowex = sum(abs(Mca[h, k]) * w[k] / w[h] for k in range(u)) / (1 - t)
            colex = sum(abs(Mca[k, h]) * w[k] / w[h] for k in range(u)) / (1 - t)
            if max(b['ROWp'], b['ROWm']) < rowex - 1e-9 or max(b['COLp'], b['COLm']) < colex - 1e-9: viol += 1
            if not (-1e-9 <= nu <= 1 + 1e-9 and -1e-9 <= theta <= 1 + 1e-9): viol += 1
print('violations', viol, ' max pi', maxpi, ' worst exact signed row', worst_ex)
# (2) grid
worst = dict(LB=np.inf, ROW=0, COL=0); arg = {}
for S in np.concatenate([np.linspace(1e-4, 1.999, 100), np.geomspace(1e-6, 1e-2, 10)]):
    for xi in np.linspace(1e-4, 1, 50):
        for theta in np.linspace(0, 1, 11):
            for nu in (0.0, 1.0):
                r = bounds(S, xi, theta, nu, one=1.0)
                row = max(r['ROWp'], r['ROWm']); col = max(r['COLp'], r['COLm'])
                if r['LB'] < worst['LB']: worst['LB'] = r['LB']; arg['LB'] = (S, xi, theta, nu)
                if row > worst['ROW']: worst['ROW'] = row; arg['ROW'] = (S, xi, theta, nu)
                if col > worst['COL']: worst['COL'] = col; arg['COL'] = (S, xi, theta, nu)
print(worst, arg)
mq = -np.inf
for xh in np.linspace(1e-4, 2, 60):
    for xk in np.linspace(1e-4, 2 - xh, 60):
        for gg in np.linspace(xh / (1 + xh) + xk / (1 + xk), 2, 30):
            mq = max(mq, pair_Q(xh, xk, gg, one=1.0))
print('max pair_Q', mq)
