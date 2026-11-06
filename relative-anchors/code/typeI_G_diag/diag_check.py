"""Check of the diagonal theorem on our own construction of the block (../typeI_G_u3/numG.py, floating point): in the
frame of numG.frame (rho_t, the merged z_t rescaled by sqrt(1 + S/sigma), the part of zeta_t orthogonal to w_U), where the
letter norm rho^2 + E|Z|^2 - |EZ|^2/(sigma + 1) is the identity, lambda_min(Ghat_tt) > (2/3)(1 - s_t/sigma) for every t,
at random weight vectors with m = 2..4 letter-2 copies and u = 3..7 free copies (spread, wide, one dominant letter-2
weight, all letter-2 weights small).  No claim rests on it; it guards the statement and its normalisation."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'typeI_G_u3'))
from numG import Ghat
rng = np.random.default_rng(11)
worst, worstcase, n = 9.0, None, 0
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
    for t in range(m):
        B = G[t * d:(t + 1) * d, t * d:(t + 1) * d]
        r = np.linalg.eigvalsh((B + B.T) / 2)[0] / (1 - sT[t] / sig)
        n += 1
        if r < worst:
            worst, worstcase = r, (m, u, round(sT[t] / sig, 4))
print(f'{n} diagonal blocks at 240 weight vectors: min lambda_min(Ghat_tt)/(1 - s_t/sigma) = {worst:.6f} at (m, u, s_t/sigma) = {worstcase}')
print('PASS' if worst > 2 / 3 else 'FAIL')
