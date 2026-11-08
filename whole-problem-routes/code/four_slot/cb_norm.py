"""Multiplicative cb-type norms of Delta_t = (Tr(X) I - t X)/(3 - t) on M_3.
||.||_cb = ||.||_diamond (Delta self-adjoint); covariant => attained at max. entangled input: ||J||_1/3.
Weyl multiplier form: Delta_t = sum_g mu_g Ad(W_g), ||Delta_t||_cb = sum|mu_g| (upper bound, equal here).
Co-cb norm ||Delta_t o T||_cb = 1 iff Choi of Delta_t o T (prop. to W_t = I - tF) is PSD iff t <= 1."""
import numpy as np
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3):
        F[3 * j + i, 3 * i + j] = 1
Om = np.eye(3).reshape(9)
for t in [1/3, 0.4, 0.5, 0.55]:
    J = (np.eye(9) - t * np.outer(Om, Om)) / (3 - t)      # Choi of Delta_t
    cb = np.abs(np.linalg.eigvalsh(J)).sum() / 3
    mu0, mug = (1/3 - t) / (3 - t), (1/3) / (3 - t)
    l1 = abs(mu0) + 8 * mug
    JT = (np.eye(9) - t * F) / (3 - t)                     # Choi of Delta_t o T
    cocb = np.abs(np.linalg.eigvalsh(JT)).sum() / 3
    print(f"t={t:.4f}: ||Delta||_cb = {cb:.6f} (Weyl l1 {l1:.6f}),  ||Delta o T||_cb = {cocb:.6f}")
