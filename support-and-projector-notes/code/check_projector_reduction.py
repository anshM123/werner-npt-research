"""Proposition 4: flattening and the two-flag projector identity, on a random rank-two C at n = 2 (so 2n = 4), t = 1/2.
  q_2n(E_1 (x) E_2 - E_2 (x) E_1) = 2 (G_11 G_22 - |G_12|^2),  G_ij = q_n(E_i, E_j),
  q_{2n+k}(P_k) = q_2n(D)/2 + 2^{-k}/4 [q(P_U) + q(P_V) + 2 (-1)^k Re q(P_U, P_V)]."""
import numpy as np
from endpoint import d, q

rng = np.random.default_rng(7)
n = 2
C = (rng.normal(size=(9, 2)) + 1j * rng.normal(size=(9, 2))) @ (rng.normal(size=(2, 9)) + 1j * rng.normal(size=(2, 9)))
U_, s_, Vh = np.linalg.svd(C)
E = [np.outer(U_[:, i], Vh[i]) for i in range(2)]          # u_i v_i^*
G = np.array([[q(E[i], n, E[j]) for j in range(2)] for i in range(2)])
D = np.kron(E[0], E[1]) - np.kron(E[1], E[0])
print(f"flattening: q_4(D) = {q(D, 2 * n).real:.12f},  2 det G = {2 * (G[0, 0] * G[1, 1] - abs(G[0, 1])**2).real:.12f}")
u, v = U_[:, :2], Vh[:2].conj().T
Uf = np.stack([np.kron(u[:, 0], u[:, 1]), np.kron(u[:, 1], u[:, 0])], 1)
Vf = np.stack([np.kron(v[:, 0], v[:, 1]), -np.kron(v[:, 1], v[:, 0])], 1)
assert np.allclose(Uf @ Vf.conj().T, D) and np.allclose(Uf.conj().T @ Uf, np.eye(2)) and np.allclose(Vf.conj().T @ Vf, np.eye(2))
print("D = U V^* with U, V isometries (flat)")
PU, PV = Uf @ Uf.conj().T, Vf @ Vf.conj().T
for k in (1, 2):
    f0 = np.zeros((d ** k, 1)); f0[0] = 1
    f1 = np.zeros((d ** k, 1)); f1[np.ravel_multi_index([1] * k, [d] * k)] = 1
    Z = (np.kron(Uf, f0) + np.kron(Vf, f1)) / np.sqrt(2)
    P = Z @ Z.conj().T
    assert np.allclose(P @ P, P) and np.isclose(np.trace(P).real, 2)
    rhs = q(D, 2 * n).real / 2 + 2.0**-k / 4 * (q(PU, 2 * n) + q(PV, 2 * n) + 2 * (-1)**k * q(PU, 2 * n, PV)).real
    print(f"k = {k}: P_k is a rank-two projector;  q(P_k) = {q(P, 2 * n + k).real:.12f},  formula = {rhs:.12f}")

# Proposition 5: Gram matrices under repeated flattening (pure algebra; the tensor-square rule was checked above)
A, B = (G[0, 0] * G[1, 1]).real, abs(G[0, 1])**2
F1, F2 = np.kron(E[0], E[1]), -np.kron(E[1], E[0])
G1 = np.array([[q(F1, 2 * n), q(F1, 2 * n, F2)], [q(F2, 2 * n, F1), q(F2, 2 * n)]])
print(f"Gram of (F_1, F_2) on 2n = {2 * n} copies: {np.round(G1.real, 12).tolist()};  predicted [[A, -B], [-B, A]] with "
      f"A = {A:.12f}, B = {B:.12f}")
