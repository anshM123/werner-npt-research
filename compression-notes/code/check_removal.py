"""Removal by compression (Section 9, "Removal by compression").
   For a rank-two C = sum_a |y_a><z_a| on n copies and a set T of copies, let S_T, S'_T be the spans of
   (I_T (x) <r|) y_a and (I_T (x) <r|) z_a over the other copies' basis vectors r.  Check the exact identity
       q_n(1/2, C) = <psi_C| (tau_T (x) rho^{(x)(n-b)})^G |psi_C>,   tau_T = (P_S (x) P_S') rho^{(x)b} (P_S (x) P_S'),
   psi_C = sum_a y_a (x) conj(z_a); and that compressing the second side to conj(S'_T) instead gives a different number.
   Then print, for random C at n = 4, the smallest normalised partial-transpose eigenvalue of tau_T for each triple."""
import itertools
import numpy as np
from werner import q

rng = np.random.default_rng(3)
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3):
        F[3 * i + j, 3 * j + i] = 1
rho = np.eye(9) - F / 2


def tens(c):
    R = np.array([[1.0]])
    for _ in range(c):
        R = np.kron(R, rho)
    perm = list(range(0, 2 * c, 2)) + list(range(1, 2 * c, 2))
    return R.reshape([3] * (4 * c)).transpose(perm + [p + 2 * c for p in perm]).reshape(9 ** c, 9 ** c)


def pt_B(M, dA, dB):
    return M.reshape(dA, dB, dA, dB).transpose(0, 3, 2, 1).reshape(dA * dB, dA * dB)


def support(Y, T, n):
    """orthonormal basis of the span of (I_T (x) <r|) y_a, as vectors on the copies of T (in increasing order)"""
    rest = [s for s in range(n) if s not in T]
    M = np.moveaxis(Y.reshape([Y.shape[0]] + [3] * n), [1 + s for s in rest], list(range(n - len(rest) + 1, n + 1)))
    M = M.reshape(Y.shape[0], 3 ** len(T), 3 ** len(rest)).transpose(0, 2, 1).reshape(-1, 3 ** len(T))
    U, s, _ = np.linalg.svd(M.T, full_matrices=False)
    return U[:, s > 1e-10 * s[0]]


cr = lambda *s: rng.standard_normal(s) + 1j * rng.standard_normal(s)
n = 4
R4 = tens(4)
Y, Z = cr(2, 81), cr(2, 81)
C = Y.T @ Z.conj()
psi = sum(np.kron(Y[a], Z[a].conj()) for a in range(2))
print(f"q_4 (dense) = {q(C, 4, 0.5):.10f};  <psi|W^G|psi> = {np.vdot(psi, pt_B(R4, 81, 81) @ psi).real:.10f}")
T = (0, 1, 2)
S, Sp = support(Y, T, n), support(Z, T, n)
for name, SB in (("P_S (x) P_S'      ", Sp), ("P_S (x) P_conj(S')", Sp.conj())):
    PA = np.kron(S @ S.conj().T, np.eye(3)); PB = np.kron(SB @ SB.conj().T, np.eye(3))
    Wc = np.kron(PA, PB) @ R4 @ np.kron(PA, PB)
    print(f"  compression {name}: <psi|(tau (x) rho)^G|psi> = {np.vdot(psi, pt_B(Wc, 81, 81) @ psi).real:.10f}")
R3 = tens(3)
print("random C at n = 4, smallest normalised PT eigenvalue of tau_T for each triple:")
for t in range(3):
    Y, Z = cr(2, 81), cr(2, 81); out = []
    for T in itertools.combinations(range(4), 3):
        S, Sp = support(Y, T, n), support(Z, T, n); W = np.kron(S, Sp); tau = W.conj().T @ R3 @ W
        kA, kB = S.shape[1], Sp.shape[1]
        out.append(np.linalg.eigvalsh(pt_B((tau + tau.conj().T) / 2, kA, kB))[0] / np.trace(tau).real * kA * kB)
    print("  " + ", ".join(f"{v:+.3f}" for v in out))
