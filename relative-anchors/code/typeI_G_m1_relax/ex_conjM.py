import numpy as np
from numlib import Gz, vvec
rng = np.random.default_rng(1)
for u in (3, 4, 6, 8):
    worst = 1e9; worstZ = -1e9; mineig = 1e9
    for it in range(40):
        st = float(np.exp(rng.uniform(-4, 4)))
        sU = np.exp(rng.uniform(-3, 3, u))
        G, C, B, A = Gz(st, sU)
        S = sU.sum()
        v = vvec(st, S, u)
        r = (G @ v) / (np.diag(G) * v)
        off = G - np.diag(np.diag(G))
        worst = min(worst, r.min()); worstZ = max(worstZ, (off/np.sqrt(np.outer(np.diag(G),np.diag(G)))).max())
        d = 1/np.sqrt(np.diag(G)); mineig = min(mineig, np.linalg.eigvalsh(G*np.outer(d,d)).min())
    print(u, 'min normalised (D) margin', worst, ' max normalised offdiag', worstZ, ' min eig of diag-normalised G', mineig)
