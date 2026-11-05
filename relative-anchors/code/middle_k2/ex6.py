"""Exploration: gap / x_min (Y units) at equal weights, k = 4..7."""
import numpy as np
from mid import Mid
from ex1 import kernel_basis
for k in (4, 5, 6, 7):
    w = np.ones(k); m = Mid(list(w), float); Y = np.array(m.full(), float) / 4
    n = m.np; N1 = k * n
    K = np.zeros((N1 + k, n + 1)); K[:N1, :n] = kernel_basis(m); K[N1:, n] = w
    Qk, _ = np.linalg.qr(K); P = np.eye(N1 + k) - Qk @ Qk.T
    ev = np.linalg.eigvalsh(P @ Y @ P + 10 * Qk @ Qk.T)
    print(f"k={k} equal weights: gap {ev[0]:.5f} (Y units), gap/x_min {ev[0]*k:.4f}, gap/ub {ev[0]/(1/(2*(k-1))):.4f}")
