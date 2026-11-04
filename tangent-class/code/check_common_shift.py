"""For general rank-two operators, the single-copy pencils can all be positive definite at one shared negative shift.

The minimiser argument of Theorem 3, run over all rank-two operators instead of the tangent class, needs at a
hypothetical negative minimiser C* that Q_s - m N_s is PSD at every copy s, with m = R(C*) < 0.  A sign-free version
would ask that no thick C has all its pencils PSD at one shared lambda in [-1/2, 0).  data/common_shift_pd_n3.npy is a
thick three-copy operator for which all three pencils are positive definite at lambda = -1/2, so that version is
false.  The coupling lambda = R(C) is what the argument actually has (here R = 0.265 > 0).
Checks: rank two; local supports of range and co-range far from thin; the 9x9 forms against direct q(GXG) for random
positive M; lambda_min(Q_s + N_s/2) > 0 at every copy; the PSD interval of every pencil."""
import numpy as np
from scipy.linalg import sqrtm
from scipy.optimize import brentq
from common import squeeze_forms, q, q_direct, Bt

rng = np.random.default_rng(0); n = 3; D = 27
X = np.load("../data/common_shift_pd_n3.npy"); X = X / np.linalg.norm(X)
S = np.linalg.svd(X, compute_uv=False)
print(f"singular values {S[0]:.4f}, {S[1]:.4f}, {S[2]:.1e}; R = {q(X, n):.5f} (direct partial traces {q_direct(X, n):.5f})")
ok = S[2] < 1e-12 * S[0] and abs(q(X, n) - q_direct(X, n)) < 1e-12


def embed(G, s):
    ops = [np.eye(3)] * n; ops[s] = G; out = ops[0]
    for o in ops[1:]:
        out = np.kron(out, o)
    return out


worst_err, worst_margin, ratios = 0, np.inf, []
for s in range(n):
    K, Kn = squeeze_forms(X, n, s)
    for _ in range(5):
        A = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); M = A @ A.conj().T; G = embed(sqrtm(M), s)
        Y = G @ X @ G; c = np.array([np.trace(B.conj().T @ M).real for B in Bt])
        worst_err = max(worst_err, abs(c @ K @ c - q_direct(Y, n)) / abs(q_direct(Y, n)),
                        abs(c @ Kn @ c - np.linalg.norm(Y) ** 2) / np.linalg.norm(Y) ** 2)
    g = lambda l: np.linalg.eigvalsh(K - l * Kn)[0]
    lo, hi = brentq(g, -100, -0.5), brentq(g, -0.5, 100)
    margin = g(-0.5) / np.linalg.norm(Kn, 2); worst_margin = min(worst_margin, margin)
    for Yp in (X @ X.conj().T, X.conj().T @ X):
        T = np.moveaxis(Yp.reshape([3] * (2 * n)), [s, n + s], [0, n]).reshape(3, D // 3, 3, D // 3)
        ev = np.linalg.eigvalsh(np.einsum('aibi->ab', T)); ratios.append(ev[0] / ev[-1])
    print(f"copy {s}: lambda_min(Q + N/2)/||N|| = {margin:+.4e}; Q - lambda N is PSD exactly for lambda in [{lo:+.4f}, {hi:+.4f}]")
print(f"forms vs direct evaluation: largest relative error {worst_err:.1e}; smallest local-support ratio (range and co-range) {min(ratios):.3f}")
ok &= worst_err < 1e-10 and worst_margin > 0 and min(ratios) > 0.25
assert ok
print("PASS")
