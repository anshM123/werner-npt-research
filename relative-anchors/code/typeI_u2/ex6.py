"""Hypothesis: the Schur complement G onto the T-blocks (rho_t, z_t, z'_t) for u >= 2 equals the u = 1 Schur complement
(with sigma = s_U) on (rho_t, z_t), and the z'_t directions decouple."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from ex5 import Gmat
np.set_printoptions(linewidth=220, precision=5, suppress=True)
rng = np.random.default_rng(2)
for u, m in [(2, 2), (3, 2), (2, 3), (4, 3)]:
    k = u + m; T = list(range(m))
    w = np.exp(rng.uniform(-1, 1, k))
    G, d, s, sT, sU, _ = Gmat(w, T, rescale=False)
    w1 = np.concatenate([w[:m], [np.sqrt(sU)]])     # u = 1 system: T = first m, U = {last}
    G1, d1, *_ = Gmat(w1, T, rescale=False)
    # pick the (rho_t, z_t) coordinates of G
    sel = [i * d + j for i in range(m) for j in range(2)]
    Gsel = G[np.ix_(sel, sel)]
    print(f"u={u} m={m}: max|G(rho,z) - G_u1| = {np.abs(Gsel - G1).max():.2e}; ")
    perp = [i * d + j for i in range(m) for j in range(2, d)]
    print("   coupling (rho,z)-(z') max", np.abs(G[np.ix_(sel, perp)]).max())
    print("   z' block:\n", G[np.ix_(perp, perp)])
