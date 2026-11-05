"""Small tools: the endpoint form q_n on dense and on sparse operators, local Weyl operators."""
import itertools
import numpy as np

OMEGA = np.exp(2j * np.pi / 3)


def q_dense(X, n):
    """q_n(X) = <vec X, M_n vec X>, M = (I - 1/2 |Omega><Omega|)^{(x) n}; X is 3^n x 3^n."""
    T = np.asarray(X, dtype=complex).reshape((3,) * n + (3,) * n)
    MT = T.copy()
    for s in range(n):
        tr = np.trace(MT, axis1=s, axis2=n + s)
        tr = np.expand_dims(np.expand_dims(tr, s), n + s)
        shape = [1] * (2 * n)
        shape[s] = 3
        shape[n + s] = 3
        MT = MT - 0.5 * np.eye(3).reshape(shape) * tr
    return np.vdot(T, MT)


def q_sparse(entries, n):
    """entries: dict {(row tuple, col tuple): value}. Uses q = sum_J (-1/2)^|J| ||Tr_J X||_F^2."""
    items = list(entries.items())
    total = 0.0
    for mask in range(1 << n):
        J = [s for s in range(n) if mask >> s & 1]
        Jc = [s for s in range(n) if not mask >> s & 1]
        acc = {}
        for (r, c), v in items:
            if all(r[s] == c[s] for s in J):
                key = (tuple(r[s] for s in Jc), tuple(c[s] for s in Jc))
                acc[key] = acc.get(key, 0) + v
        total += (-0.5) ** len(J) * sum(abs(v) ** 2 for v in acc.values())
    return total


def weyl1(a, b):
    """X^a Z^b on one qutrit: |j> -> omega^{b j} |j + a>."""
    W = np.zeros((3, 3), dtype=complex)
    for j in range(3):
        W[(j + a) % 3, j] = OMEGA ** (b * j)
    return W


def weyl(labels):
    W = np.array([[1.0 + 0j]])
    for (a, b) in labels:
        W = np.kron(W, weyl1(a, b))
    return W


def local_weyl(k, lam):
    """U_lambda = (X^a Z^b)^{(x) k}."""
    return weyl([lam] * k)


def rand_c(shape, rng):
    return rng.standard_normal(shape) + 1j * rng.standard_normal(shape)


def marginals_sparse(vec):
    """vec: dict {string tuple: amplitude}. Returns the list of one-copy marginals (3 x 3)."""
    n = len(next(iter(vec)))
    out = []
    for s in range(n):
        groups = {}
        for x, v in vec.items():
            groups.setdefault(x[:s] + x[s + 1:], []).append((x[s], v))
        rho = np.zeros((3, 3), dtype=complex)
        for lst in groups.values():
            for (i, u) in lst:
                for (j, w) in lst:
                    rho[i, j] += u * np.conj(w)
        out.append(rho)
    return out


def q_sparse_np(R, Cc, V):
    """Vectorised q_n for a sparse operator: R, Cc integer arrays (nnz x n) of row/column letters, V values."""
    R = np.asarray(R); Cc = np.asarray(Cc); V = np.asarray(V, dtype=complex)
    n = R.shape[1]
    total = 0.0
    eq = (R == Cc)
    for mask in range(1 << n):
        J = [s for s in range(n) if mask >> s & 1]
        Jc = [s for s in range(n) if not mask >> s & 1]
        sel = np.all(eq[:, J], axis=1) if J else np.ones(len(V), dtype=bool)
        if not sel.any():
            continue
        if Jc:
            w = 3 ** np.arange(len(Jc))
            key = (R[sel][:, Jc] @ w) * (3 ** len(Jc)) + (Cc[sel][:, Jc] @ w)
        else:
            key = np.zeros(sel.sum(), dtype=np.int64)
        _, inv = np.unique(key, return_inverse=True)
        sums = np.zeros(inv.max() + 1, dtype=complex)
        np.add.at(sums, inv, V[sel])
        total += (-0.5) ** len(J) * float(np.sum(np.abs(sums) ** 2))
    return total


def dense_to_sparse(X, n):
    idx = np.argwhere(np.abs(X) > 0)
    R = np.array([[int(i) // 3 ** (n - 1 - s) % 3 for s in range(n)] for i in idx[:, 0]])
    Cc = np.array([[int(j) // 3 ** (n - 1 - s) % 3 for s in range(n)] for j in idx[:, 1]])
    return R, Cc, X[idx[:, 0], idx[:, 1]]
