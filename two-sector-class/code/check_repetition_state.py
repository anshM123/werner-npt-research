"""LOCC reading of the three-sector collapse.  Project m (odd) copies of the endpoint Werner state locally onto the
repetition code span{1..1, 2..2, 3..3} on both sides.  The partial transpose of the result is (up to normalisation)
2^{1-m} Z_1 + (1 - 2^{1-m}) D,  Z_1 = 1 - |Omega><Omega|/2 (one copy), D = sum_{i != j} |ij><ij|.
Checks: the compression identity, positivity of the state, its PT spectrum (smallest eigenvalue -2^{-m}), and the
parameter of its Werner twirl, a = -(1 + 2^{1-m})/(3 + 2^{1-m}).   usage: python check_repetition_state.py m"""
import sys, numpy as np
m = int(sys.argv[1]); al = 2.0 ** (1 - m)
Om = np.eye(3).reshape(9); Z1 = np.eye(9) - 0.5 * np.outer(Om, Om)
Zm = np.ones((1, 1))
for _ in range(m): Zm = np.kron(Zm, Z1)
Zm = Zm.reshape([3, 3] * m + [3, 3] * m)
perm = [2 * s for s in range(m)] + [2 * s + 1 for s in range(m)]; perm = perm + [2 * m + p for p in perm]
Zm = Zm.transpose(perm).reshape(3 ** (2 * m), 3 ** (2 * m))
G = np.zeros((3 ** m, 3))
for i in range(3): G[int(np.ravel_multi_index((i,) * m, [3] * m)), i] = 1
V = np.kron(G, G); Zlog = V.T @ Zm @ V
D = np.diag([0. if i == j else 1. for i in range(3) for j in range(3)])
pred = al * Z1 + (1 - al) * D
rho = pred.reshape(3, 3, 3, 3).transpose(0, 3, 2, 1).reshape(9, 9)
F = np.eye(9).reshape(3, 3, 3, 3).transpose(0, 1, 3, 2).reshape(9, 9)
a_tw = np.trace(rho @ F) / np.trace(rho)            # Werner twirl: rho -> (1 + a F)/(9 + 3a) with Tr(rho F)/Tr(rho) = (3 + 9a)/(9 + 3a)
a = (9 * a_tw - 3) / (9 - 3 * a_tw)
print(f"m={m}: compression error {np.abs(Zlog - pred).max():.1e}; state min eigenvalue {np.linalg.eigvalsh(rho).min():+.4f}; "
      f"PT min eigenvalue {np.linalg.eigvalsh(pred).min():+.6f} (-2^-m = {-2.0 ** -m:+.6f}); Werner twirl a = {a:+.6f} "
      f"(formula {-(1 + al) / (3 + al):+.6f})")
