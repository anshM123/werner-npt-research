"""The sufficient condition left by Corollary G' (Vac(G) + ||Gamma_2 G_2||^2 >= ||Gamma_1 G_1||^2 for all G) fails:
lambda_max of Gamma_1^* Gamma_1 - Vac_11 on H_{n-1} (G_2 = 0) at random tangent data, basis a = |a| e_1."""
import numpy as np
from fock_tangent import *
rng = np.random.default_rng(9); worst = -1e9
for n in (3, 4):
    for trial in range(5):
        U, a = rand_tangent(rng, n)
        W = np.array([[a[0], -np.conj(a[1])], [a[1], np.conj(a[0])]]) / np.linalg.norm(a); U = [u @ np.conj(W) for u in U]; a = W.conj().T @ a
        O = Ops(U, np.array([abs(a[0]), 0])); D = O.D
        idx = O.sector_idx(n - 1)[:len(O.sector_idx(n - 1)) // 2]            # reference-1 copies of sector n-1 strings
        sel = idx                                                            # G_1 lives at r = 1 (first block)
        Gam1sq = (O.al.conj().T)[:D, :][:, :]                                # alpha^* : strings -> C^2 (x) strings; r=1 block
        M = (O.al.conj().T)[:D, sel]                                         # Gamma_1 restricted to sector n-1 strings
        G1 = M.conj().T @ M
        Vac11 = O.Vac[np.ix_(sel, sel)]
        lm = np.linalg.eigvalsh(G1 - Vac11)[-1]; worst = max(worst, lm)
        print(f"n={n} trial {trial}: lambda_max(Gamma_1^*Gamma_1 - Vac_11) = {lm:+.4f}")
print("largest", worst); assert worst > 1e-3; print("PASS (the sufficient condition fails, as claimed)")
