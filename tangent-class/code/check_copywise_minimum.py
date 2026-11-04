"""Tangent operators that are squeeze minima at every copy (against Conditions SC-cw, SC-desc and SC'-infinity).

data/copywise_squeeze_min_n3.npy and data/copywise_squeeze_min_n4.npy have range exactly inside T_e (e = |0...0>), rank
two, and range and co-range local-support ratios >= 0.05 at every copy.  At every copy s, F_s = Q_s - R(C) N_s is
positive semidefinite, its only null direction is H = I, and F_s(I, .) = 0 (C is critical).  So C minimises R over
the squeeze orbit of each copy; this is also checked directly on sampled squeezes.  The pencils are PSD only for
lambda in a short interval just below R(C) > 0, so these points do not violate Condition SC'."""
import numpy as np, itertools
from scipy.linalg import sqrtm, expm
from scipy.optimize import brentq
from common import squeeze_forms, q, q_direct, Bt

rng = np.random.default_rng(1); ok = True
cI = np.array([np.trace(B.conj().T @ np.eye(3)).real for B in Bt]); cI /= np.linalg.norm(cI)
for name in ("copywise_squeeze_min_n3", "copywise_squeeze_min_n4"):
    X = np.load(f"../data/{name}.npy"); n = int(round(np.log(X.shape[0]) / np.log(3))); D = 3 ** n
    X /= np.linalg.norm(X); R0 = q(X, n)
    st = list(itertools.product(range(3), repeat=n)); outside = [i for i, s in enumerate(st) if sum(c != 0 for c in s) > 1]
    Uu, S, _ = np.linalg.svd(X); leak = np.linalg.norm(Uu[outside, :2])
    print(f"{name}: R = {R0:.5f} (direct {q_direct(X, n):.5f}); rank {np.linalg.matrix_rank(X, 1e-9)}; range weight outside T_e {leak:.1e}")
    ok &= leak < 1e-10 and np.linalg.matrix_rank(X, 1e-9) == 2
    for s in range(n):
        K, Kn = squeeze_forms(X, n, s); F = K - R0 * Kn; w, V = np.linalg.eigh(F); sc = np.linalg.norm(Kn, 2)
        crit = np.linalg.norm(F @ cI) / sc; nullI = abs(V[:, 0] @ cI)
        rat = []
        for Y in (X @ X.conj().T, X.conj().T @ X):
            T = np.moveaxis(Y.reshape([3] * (2 * n)), [s, n + s], [0, n]).reshape(3, D // 3, 3, D // 3)
            e = np.linalg.eigvalsh(np.einsum('aibi->ab', T)); rat.append(e[0] / e[-1])
        def emb(G):
            ops = [G if t == s else np.eye(3) for t in range(n)]; out = ops[0]
            for o in ops[1:]:
                out = np.kron(out, o)
            return out
        worst = np.inf
        for k in range(40):
            A = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); M = expm((0.05 if k < 20 else 1.0) * (A + A.conj().T) / 2)
            G = emb(sqrtm(M)); Y = G @ X @ G; worst = min(worst, q(Y / np.linalg.norm(Y), n) - R0)
        g = lambda l: np.linalg.eigvalsh(K - l * Kn)[0] / sc
        lo = brentq(g, -50, R0 - 1e-9)
        print(f"   copy {s}: ratios {rat[0]:.3f}/{rat[1]:.3f}; eigenvalues of F/||N||: min {w[0] / sc:+.1e} (null vector . I = {nullI:.6f}), "
              f"next {w[1] / sc:+.2e}; |F(I, .)| {crit:.1e}; min sampled R(GCG) - R(C) {worst:+.1e}; pencil PSD on [{lo:.4f}, {R0:.4f}]")
        ok &= abs(w[0] / sc) < 1e-10 and w[1] / sc > 1e-5 and nullI > 0.9999 and crit < 1e-8 and worst > -1e-10 and min(rat) > 0.045 and lo > 0
assert ok
print("PASS")
