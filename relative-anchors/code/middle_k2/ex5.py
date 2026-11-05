"""Exploration: as one weight -> 0, gap / ub -> 1 and the soft eigenvectors -> Psi(a,h) (overlap)."""
import numpy as np
from mid import Mid
from ex1 import kernel_basis
from ex2 import soft_basis
for k in (4, 5):
    for eps in (0.3, 0.1, 0.03, 0.01):
        w = np.array([1.0, 1.3, 0.8, 1.1, 0.9][:k]); w[-1] = eps
        m = Mid(list(w), float); Y = np.array(m.full(), float) / 4
        n = m.np; N1 = k * n; s = w**2; S = s.sum()
        K = np.zeros((N1 + k, n + 1)); K[:N1, :n] = kernel_basis(m); K[N1:, n] = w
        Qk, _ = np.linalg.qr(K); P = np.eye(N1 + k) - Qk @ Qk.T
        ev, V = np.linalg.eigh(P @ Y @ P + 10 * Qk @ Qk.T)
        ub = s[-1] / (2 * (S - s[-1]))
        Sh = np.linalg.qr(soft_basis(m)[k - 1])[0]
        ov = [np.linalg.norm(Sh.T @ V[:N1, i]) ** 2 for i in range(k - 2)]
        print(f"k={k} eps={eps}: lowest {np.round(ev[:k-1], 6)}, ub {ub:.6f}, ratio {np.round(ev[:k-2]/ub, 4)}, overlap with S_h {np.round(ov, 4)}")
