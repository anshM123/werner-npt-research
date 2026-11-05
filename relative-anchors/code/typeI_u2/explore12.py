import numpy as np
from rowbounds import bounds
worst = dict(LBc=np.inf, ROW=0, COL=0)
arg = {}
for S in np.concatenate([np.linspace(1e-4, 1.999, 120), np.geomspace(1e-6, 1e-2, 20)]):
    for xi in np.linspace(1e-4, 1, 60):
        for theta in np.linspace(0, 1, 15):
            for nu in (0.0, 1.0):
                r = bounds(S, xi, theta, nu, one=1.0)
                lbc = r['LBc_a'] if r['mh'] >= 0 else r['LBc_b']
                row = max(r['ROW_p'], r['ROW_m']); col = max(r['COL_p'], r['COL_m'])
                if lbc < worst['LBc']: worst['LBc'] = lbc; arg['LBc'] = (S, xi, theta, nu)
                if row > worst['ROW']: worst['ROW'] = row; arg['ROW'] = (S, xi, theta, nu)
                if col > worst['COL']: worst['COL'] = col; arg['COL'] = (S, xi, theta, nu)
print(worst); print(arg)
