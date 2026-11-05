"""Exploration (floating point; no claim rests on it): equal U-weights, how Ghat depends on u.

At s_h = S/u for all h in U the symmetric group of U acts, so Ghat splits into the trivial part (rho_t, c z_t) and
u - 1 copies of the standard part (zeta'_t).  We print, for fixed T-weights and fixed S, the trivial 2m x 2m block and
the standard m x m block for u = 1..12 (u = 1 is the merged system), and the distance between consecutive u."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from numG import Ghat

np.set_printoptions(precision=5, suppress=True, linewidth=160)
for sT, S in (([0.3, 0.5, 0.2], 1.0), ([0.05, 0.02, 0.01], 1.0), ([2.0, 1.0], 1.0), ([0.004], 1.0)):
    m = len(sT)
    print(f'sT = {sT}, S = {S}')
    prev = None
    for u in list(range(1, 9)) + [12]:
        sU = [S / u] * u
        Gh = Ghat(sT, sU)
        d = u + 1
        mz = [t * d + i for t in range(m) for i in (0, 1)]
        A = Gh[np.ix_(mz, mz)]
        line = f'  u={u:2d}: lam_min(triv) = {np.linalg.eigvalsh(A).min():.5f}'
        if u >= 2:
            perp = [t * d + 2 for t in range(m)]
            B = Gh[np.ix_(perp, perp)]
            line += f', std block eig = {np.round(np.linalg.eigvalsh(B), 5)}, coupling = {np.abs(Gh[np.ix_(mz, perp)]).max():.1e}'
        if prev is not None:
            line += f', |triv(u) - triv(u_prev)| = {np.abs(A - prev).max():.5f}'
        prev = A
        print(line)
        if u in (1, 2, 12):
            print('     triv block:\n', A)
