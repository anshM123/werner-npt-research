"""Tangent operators with every one-copy pencil positive definite at one common negative lambda (against SC').

data/shared_shift_pd_T_n3.npy and data/shared_shift_pd_T_n4.npy have range exactly inside T_e (e = |0...0>), rank two,
and range and co-range local-support ratios >= 0.05 at every copy.  At lambda = -0.01, Q_s - lambda N_s is positive
definite at every copy s.  The forms are checked against direct evaluation of q(GCG) at random positive M."""
import numpy as np, itertools
from scipy.linalg import sqrtm
from common import squeeze_forms, q, q_direct, Bt

rng = np.random.default_rng(3); ok = True; lam = -0.01
for name in ("shared_shift_pd_T_n3", "shared_shift_pd_T_n4"):
    X = np.load(f"../data/{name}.npy"); n = int(round(np.log(X.shape[0]) / np.log(3))); D = 3 ** n; X /= np.linalg.norm(X)
    st = list(itertools.product(range(3), repeat=n)); outside = [i for i, s in enumerate(st) if sum(c != 0 for c in s) > 1]
    Uu, S, _ = np.linalg.svd(X); leak = np.linalg.norm(Uu[outside, :2]); rk = np.linalg.matrix_rank(X, 1e-9)
    print(f"{name}: R = {q(X, n):.5f}; rank {rk}; range weight outside T_e {leak:.1e}")
    ok &= leak < 1e-10 and rk == 2
    for s in range(n):
        K, Kn = squeeze_forms(X, n, s); sc = np.linalg.norm(Kn, 2); m = np.linalg.eigvalsh(K - lam * Kn)[0] / sc
        lams = np.linspace(-0.6, 0.6, 2401); pd = [l for l in lams if np.linalg.eigvalsh(K - l * Kn)[0] > 0]
        A = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); M = A @ A.conj().T
        ops = [sqrtm(M) if t == s else np.eye(3) for t in range(n)]; G = ops[0]
        for o in ops[1:]:
            G = np.kron(G, o)
        Y = G @ X @ G; c = np.array([np.trace(B.conj().T @ M).real for B in Bt]); err = abs(c @ K @ c - q_direct(Y, n)) / abs(q_direct(Y, n))
        rat = []
        for Yp in (X @ X.conj().T, X.conj().T @ X):
            T = np.moveaxis(Yp.reshape([3] * (2 * n)), [s, n + s], [0, n]).reshape(3, D // 3, 3, D // 3)
            e = np.linalg.eigvalsh(np.einsum('aibi->ab', T)); rat.append(e[0] / e[-1])
        print(f"   copy {s}: lambda_min(Q - lambda N)/||N|| at lambda = {lam} is {m:+.3e}; PD for lambda in about [{min(pd):+.3f}, {max(pd):+.3f}]; "
              f"ratios {rat[0]:.3f}/{rat[1]:.3f}; forms vs direct {err:.1e}")
        ok &= m > 1e-5 and err < 1e-10 and min(rat) > 0.045
assert ok
print("PASS")
