"""A tangent operator that is a strict minimum of R on the squeeze orbit of one copy.

data/percopy_squeeze_min_n3.npy is in the tangent class (found by code/search_strong.py, mode min0), has rank two and
is thick at every copy (range and co-range local-support ratios >= 0.15).  At copy 0, F = Q_0 - R(C) N_0 is positive
semidefinite: its only null direction is H = I (the trivial rescaling), and the other eight eigenvalues are positive.
So C is a minimum of R over all squeezes of copy 0, and the per-copy form of the one-copy condition is false.  At
copies 1 and 2, F has negative eigenvalues, so C is not a minimum there (Conditions SC-cw and SC-desc hold here)."""
import numpy as np
from common import squeeze_forms, q, Bt

n = 3; D = 27
X = np.load("../data/percopy_squeeze_min_n3.npy"); X /= np.linalg.norm(X); R0 = q(X, n)
print(f"R = {R0:.5f}; singular values {np.round(np.linalg.svd(X, compute_uv=False)[:3], 4)}")
ok = np.linalg.matrix_rank(X, 1e-10) == 2
cI = np.array([np.trace(B.conj().T @ np.eye(3)).real for B in Bt]); cI /= np.linalg.norm(cI)
for s in range(n):
    K, Kn = squeeze_forms(X, n, s); F = K - R0 * Kn; w, V = np.linalg.eigh(F); w = w / np.linalg.norm(Kn, 2)
    rat = []
    for Y in (X @ X.conj().T, X.conj().T @ X):
        T = np.moveaxis(Y.reshape([3] * 6), [s, 3 + s], [0, 3]).reshape(3, 9, 3, 9); e = np.linalg.eigvalsh(np.einsum('aibi->ab', T)); rat.append(e[0] / e[-1])
    print(f"copy {s}: local-support ratios {rat[0]:.3f} (range), {rat[1]:.3f} (co-range); eigenvalues of (Q - R N)/||N||: "
          f"{np.array2string(w, precision=4, max_line_width=250)}")
    ok &= min(rat) > 0.14
    if s == 0:
        overlap = abs(V[:, 0] @ cI)
        print(f"   copy 0: null direction overlaps H = I with |cos| = {overlap:.6f}; second eigenvalue {w[1]:+.3e}")
        ok &= abs(w[0]) < 1e-8 and w[1] > 1e-5 and overlap > 0.999
    else:
        ok &= w[0] < -1e-3
assert ok
print("PASS")
