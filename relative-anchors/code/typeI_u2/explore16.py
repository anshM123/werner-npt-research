import numpy as np
from rowbounds2 import bounds
for tau in (0.72, 0.75, 0.8):
    worstI = np.inf; worstII = 0
    for S in np.concatenate([np.linspace(1e-4, 1.9999, 160), np.geomspace(1e-7, 1e-2, 12)]):
        for xi in np.concatenate([np.linspace(1e-5, 1, 60), [1e-7]]):
            for theta in np.linspace(0, 1, 11):
                for nu in (0.0, 1.0):
                    r = bounds(S, xi, theta, nu, one=1.0)
                    v = r['LB'] - tau * max(r['ROWp'], r['ROWm'])
                    if v < worstI: worstI = v; argI = (S, xi, theta, nu)
                    c = max(r['COLp'], r['COLm'])
                    if c > worstII: worstII = c; aII = (S, xi, theta, nu)
    print('tau', tau, ' min LB - tau ROW', round(worstI, 4), argI, '  max COL', round(worstII, 4), aII)
