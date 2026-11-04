"""Covariance K^{UM}(F) = K^U(conj(M) F) of the single-excitation form (used to pass from normalised to arbitrary U)."""
import numpy as np
from fock_tangent import *
rng = np.random.default_rng(5); worst = 0
for n in (3, 4):
    for t in range(3):
        U = [rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)) for _ in range(n)]; M = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        O1 = Ops(U, np.zeros(2)); O2 = Ops([u @ M for u in U], np.zeros(2)); D = O1.D
        K = lambda O: O.AbU + O.Vac - O.T + 0.5 * O.al.conj().T @ O.al
        F = rng.normal(size=2 * D) + 1j * rng.normal(size=2 * D); G = np.kron(np.conj(M), np.eye(D)) @ F
        v1 = np.vdot(F, K(O2) @ F).real; v2 = np.vdot(G, K(O1) @ G).real; worst = max(worst, abs(v1 - v2) / abs(v1))
print("covariance worst rel err", worst); assert worst < 1e-10; print("PASS")
