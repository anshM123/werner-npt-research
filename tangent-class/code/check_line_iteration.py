"""Sanity check of Theorem 3' on examples: at each thick copy take a Finsler direction H (N_s(H) = 0, Q_s(H) < 0, from a
small SDP), walk the line M = I + eps H to both ends of its PSD interval, and recurse on the end operators until every
copy is thin.  Checks: R is quasi-concave on each line (R(0) >= min of the end values), and all leaves have q >= 0.
Needs cvxpy (solver CLARABEL)."""
import numpy as np, cvxpy as cp
from scipy.linalg import sqrtm
from common import tangent_range, squeeze_forms, q, Bt, d
rng = np.random.default_rng(5)


def local_dims(X, n, tol=1e-7):
    Uu, S, _ = np.linalg.svd(X); r = (S > tol * S[0]).sum(); u = (Uu[:, :r] * S[:r]).T
    R = sum(np.outer(u[k], u[k].conj()) for k in range(r)).reshape([d] * (2 * n)); dims = []
    for s in range(n):
        T = np.moveaxis(R, [s, n + s], [0, n]).reshape(d, d ** (n - 1), d, d ** (n - 1)); ev = np.linalg.eigvalsh(np.einsum('aibi->ab', T))
        dims.append(int((ev > tol * ev[-1]).sum()))
    return dims


def finsler(K, Kn):
    Z = cp.Variable((9, 9), symmetric=True)
    pr = cp.Problem(cp.Minimize(cp.trace(Z @ K)), [Z >> 0, cp.trace(Z @ Kn) == 0, cp.trace(Z) == 1]); pr.solve(solver=cp.CLARABEL)
    w, V = np.linalg.eigh(Z.value); h = V[:, -1]
    # top eigenvector of the SDP solution; Z is rank one up to solver accuracy, so N(H) = 0 up to that accuracy
    return np.einsum('k,kab->ab', h, Bt), pr.value


def squeeze(X, n, s, M):
    G = sqrtm(M); G = (G + G.conj().T) / 2
    T = X.reshape([d] * (2 * n)); T = np.moveaxis(np.tensordot(G, T, axes=([1], [s])), 0, s)
    T = np.moveaxis(np.tensordot(T, G, axes=([n + s], [0])), -1, n + s); return T.reshape(d ** n, d ** n)


stats = {"lines": 0, "worst_drop": np.inf, "leaves": 0, "min_leaf": np.inf}


def descend(X, n, depth=0):
    X = X / np.linalg.norm(X); dims = local_dims(X, n); R0 = q(X, n)
    thick = [s for s in range(n) if dims[s] == 3]
    if not thick or depth > n:
        stats["leaves"] += 1; stats["min_leaf"] = min(stats["min_leaf"], R0); return R0
    s = thick[0]; K, Kn = squeeze_forms(X, n, s); H, val = finsler(K, Kn)
    ev = np.linalg.eigvalsh(H); em, ep = -1 / ev[-1], -1 / ev[0]                 # PSD interval of I + eps H
    ends = []
    for e in (em * (1 - 1e-9), ep * (1 - 1e-9)):
        Y = squeeze(X, n, s, np.eye(3) + e * H)
        ends.append(descend(Y, n, depth + 1) if np.linalg.norm(Y) > 1e-12 else np.inf)
    stats["lines"] += 1; stats["worst_drop"] = min(stats["worst_drop"], R0 - min(ends))
    return min(ends)


for n, trials in ((3, 3), (4, 1)):
    D = d ** n
    for trial in range(trials):
        U = [rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)) for _ in range(n)]; a = rng.normal(size=2) + 1j * rng.normal(size=2)
        u = tangent_range(U, a, n); v = rng.normal(size=(2, D)) + 1j * rng.normal(size=(2, D))
        X = np.outer(u[0], v[0].conj()) + np.outer(u[1], v[1].conj())
        R0 = q(X / np.linalg.norm(X), n); m = descend(X, n)
        print(f"n={n} trial {trial}: R(C) = {R0:.4f}, minimum over the leaves = {m:.4f}")
print(f"lines followed {stats['lines']}, smallest R(start) - min(R(ends)) {stats['worst_drop']:+.2e} (>= 0: quasi-concavity), "
      f"leaves {stats['leaves']}, smallest leaf value {stats['min_leaf']:+.4f}")
assert stats["worst_drop"] > -1e-7 and stats["min_leaf"] >= 0
print("PASS")
