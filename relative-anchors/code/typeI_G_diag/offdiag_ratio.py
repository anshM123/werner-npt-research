"""Floating-point look at the open half of the weighted Gershgorin route, on our own block (../typeI_G_u3/numG.py): with
the diagonal theorem, G > 0 would follow from ||Ghat_ij||_op <= (3/4) sqrt(s_i s_j)/sigma for i != j.  This prints the
largest ratio ||Ghat_ij||_op / (sqrt(s_i s_j)/sigma) at random weight vectors.  It proves nothing; it shows whether the
route can close with the constant 3/4."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'typeI_G_u3'))
from numG import Ghat
rng = np.random.default_rng(12)
worst, wc, n = 0.0, None, 0
for it in range(240):
    m = int(rng.integers(2, 5)); u = int(rng.integers(3, 8)); mode = it % 4
    if mode == 0:
        sT = list(np.exp(rng.normal(0, 1.5, m))); sU = list(np.exp(rng.normal(0, 1.5, u)))
    elif mode == 1:
        sT = list(np.exp(rng.uniform(-6, 6, m))); sU = list(np.exp(rng.uniform(-6, 6, u)))
    elif mode == 2:
        sT = list(np.exp(rng.normal(0, 1, m))); sT[0] *= np.exp(rng.uniform(4, 8)); sU = list(np.exp(rng.normal(0, 1, u)))
    else:
        sT = list(np.exp(rng.normal(-5, 1, m))); sU = list(np.exp(rng.normal(0, 2, u)))
    G = Ghat(sT, sU); d = u + 1; sig = sum(sT)
    for i in range(m):
        for j in range(i + 1, m):
            B = G[i * d:(i + 1) * d, j * d:(j + 1) * d]
            r = np.linalg.norm(B, 2) / (np.sqrt(sT[i] * sT[j]) / sig)
            n += 1
            if r > worst:
                worst, wc = r, (m, u, round(sT[i] / sig, 4), round(sT[j] / sig, 4))
print(f'{n} off-diagonal blocks: max ||Ghat_ij||_op / (sqrt(s_i s_j)/sigma) = {worst:.6f} at (m, u, s_i/sigma, s_j/sigma) = {wc}')
print('the route with constant 3/4 would need this to stay below 3/4')
