"""Gelfand-Tsetlin (Young orthogonal form) bases for S_k acting on (R^m)^{(x)k} by permuting tensor factors.

gt_bases(m, k) returns {lam: {'tabs': [T...], 'W': {T: (m^k, d) array}, 'dsp': #tableaux, 'd': d}} where, for every
standard tableau T of shape lam (at most m rows), the columns of W[T] form an orthonormal basis of V_lam (x) e_T, with the
SAME basis of V_lam for every T (consistent through Young's orthogonal form).  Then for any operator X,
    sum_T W[T]^T X W[T] = dsp * (the lam-block of the S_k-twirl of X),
and an S_k-invariant Z is PSD iff each lam-block is.  Self-test at the bottom."""
import itertools
import numpy as np


def partitions(k, maxrows):
    out = []

    def rec(rem, maxpart, cur):
        if rem == 0:
            out.append(tuple(cur))
            return
        if len(cur) == maxrows:
            return
        for p in range(min(rem, maxpart), 0, -1):
            rec(rem - p, p, cur + [p])
    rec(k, k, [])
    return out


def tableaux(lam):
    """standard Young tableaux of shape lam, as tuples pos[i] = (row, col) of entry i (0-based entries)"""
    k = sum(lam)
    res = []

    def rec(i, rows, pos):
        if i == k:
            res.append(tuple(pos))
            return
        for r in range(len(lam)):
            c = rows[r]
            if c < lam[r] and (r == 0 or rows[r - 1] > c):
                rows[r] += 1
                rec(i + 1, rows, pos + [(r, c)])
                rows[r] -= 1
    rec(0, [0] * len(lam), [])
    return res


def contents(T):
    return [c - r for (r, c) in T]


def perm_matrix(m, k, perm):
    """matrix of the permutation of tensor factors: (P v)[i_{perm(0)} ... ] -- factor s goes to position perm[s]"""
    D = m ** k
    P = np.zeros((D, D))
    for idx in itertools.product(range(m), repeat=k):
        new = [0] * k
        for s in range(k):
            new[perm[s]] = idx[s]
        a = 0
        b = 0
        for s in range(k):
            a = a * m + idx[s]
            b = b * m + new[s]
        P[b, a] = 1
    return P


def transposition(m, k, i, j):
    perm = list(range(k))
    perm[i], perm[j] = j, i
    return perm_matrix(m, k, perm)


def gt_bases(m, k, seed=0):
    if k == 0:
        return {(): {'tabs': [()], 'W': {(): np.ones((1, 1))}, 'dsp': 1, 'd': 1}}
    D = m ** k
    rng = np.random.default_rng(seed)
    trans = {}
    for i in range(k):
        for j in range(i + 1, k):
            trans[(i, j)] = transposition(m, k, i, j)
    JM = [sum((trans[(j, i)] for j in range(i)), np.zeros((D, D))) for i in range(k)]   # X_0 = 0
    cvec = rng.standard_normal(k)
    Y = sum(c * X for c, X in zip(cvec, JM))
    ev, evec = np.linalg.eigh((Y + Y.T) / 2)
    out = {}
    for lam in partitions(k, m):
        tabs = tableaux(lam)
        W = {}
        T0 = tabs[0]
        target = float(np.dot(cvec, contents(T0)))
        sel = np.abs(ev - target) < 1e-7
        W[T0] = evec[:, sel]
        d = W[T0].shape[1]
        # BFS via adjacent transpositions s_i = (i, i+1)
        frontier = [T0]
        while frontier:
            new = []
            for T in frontier:
                cT = contents(T)
                for i in range(k - 1):
                    T2 = list(T)
                    T2[i], T2[i + 1] = T2[i + 1], T2[i]
                    T2 = tuple(T2)
                    if T2 in W or T2 not in tabs:
                        continue
                    r = cT[i + 1] - cT[i]          # axial distance
                    s = trans[(i, i + 1)]
                    W[T2] = (s @ W[T] - (1.0 / r) * W[T]) / np.sqrt(1 - 1.0 / r ** 2)
                    new.append(T2)
            frontier = new
        assert len(W) == len(tabs), (lam, len(W), len(tabs))
        out[lam] = {'tabs': tabs, 'W': W, 'dsp': len(tabs), 'd': d}
    return out


if __name__ == '__main__':
    for m, k in [(2, 3), (3, 3), (2, 4), (3, 4), (3, 5), (2, 6), (3, 6)]:
        B = gt_bases(m, k)
        Wall = np.hstack([B[l]['W'][T] for l in B for T in B[l]['tabs']])
        D = m ** k
        assert Wall.shape == (D, D), Wall.shape
        err_orth = np.abs(Wall.T @ Wall - np.eye(D)).max()
        # an S_k-invariant operator: A^{(x)k} with random symmetric A; blocks must be T-independent and off-T zero
        A = np.random.default_rng(1).standard_normal((m, m))
        A = A + A.T
        X = np.array([[1.0]])
        for _ in range(k):
            X = np.kron(X, A)
        err_inv = 0.0
        for l in B:
            blocks = [B[l]['W'][T].T @ X @ B[l]['W'][T] for T in B[l]['tabs']]
            err_inv = max(err_inv, max(np.abs(b - blocks[0]).max() for b in blocks))
            T0, T1 = B[l]['tabs'][0], B[l]['tabs'][-1]
            if T0 != T1:
                err_inv = max(err_inv, np.abs(B[l]['W'][T0].T @ X @ B[l]['W'][T1]).max())
        print(f'm = {m}, k = {k}: irreps {[(l, B[l]["d"], B[l]["dsp"]) for l in B]}  orth err {err_orth:.1e}  '
              f'invariance err {err_inv:.1e}', flush=True)
