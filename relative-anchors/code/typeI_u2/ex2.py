"""Exploration: comparison-matrix test of block dominance after (1) eliminating X_T exactly, (2) the gauge
w_U^T X_U = 0, split X_U = X' + a w_U^T (U-hard, U-soft), (3) optional rescaling of the soft directions."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cf import TypeI
np.set_printoptions(linewidth=200, precision=4, suppress=True)


def null_basis(A, tol=1e-12):
    u, s, vt = np.linalg.svd(A)
    r = (s > tol * max(1, s.max())).sum()
    return vt[r:].T


def reduced(w, T):
    ti = TypeI(w, T, F=float)
    M = np.array(ti.build())
    k = len(w); U = ti.U; u = len(U); w = np.asarray(w, float)
    iXT = [ti.idx[('P', t, h)] for t in T for h in U]
    rest = [i for i in range(ti.n) if i not in iXT]
    S = M[np.ix_(rest, rest)] - M[np.ix_(rest, iXT)] @ np.linalg.solve(M[np.ix_(iXT, iXT)], M[np.ix_(iXT, rest)])
    labs = [ti.labels[i] for i in rest]
    pos = {l: i for i, l in enumerate(labs)}
    wU = w[U]
    # coordinates: X_U as u x u (row r, col h)
    def XUvec(X):
        v = np.zeros(len(labs))
        for a, r in enumerate(U):
            for b, h in enumerate(U):
                v[pos[('P', r, h)]] = X[a, b]
        return v
    # U-hard basis: X' with X' w = 0 and w^T X' = 0
    Bw = null_basis(wU[None, :])            # u x (u-1), orthonormal, perp to w
    hard = [XUvec(np.outer(Bw[:, i], Bw[:, j])) for i in range(u - 1) for j in range(u - 1)]
    soft = [XUvec(np.outer(Bw[:, i], wU / np.linalg.norm(wU))) for i in range(u - 1)]
    qv = np.zeros(len(labs)); qv[pos[('Q',)]] = 1
    blocks = {'Uh': hard + [qv], 'Us': soft}
    for t in T:
        rv = np.zeros(len(labs)); rv[pos[('R', t)]] = 1
        zs = []
        for i in range(u - 1):
            z = np.zeros(len(labs))
            for b, h in enumerate(U):
                z[pos[('Z', t, h)]] = Bw[b, i]
            zs.append(z)
        zw = np.zeros(len(labs))
        for b, h in enumerate(U):
            zw[pos[('Z', t, h)]] = wU[b] / np.linalg.norm(wU)
        blocks[('T', t)] = [rv, zw] + zs
    # check kernel: kernel of S should be X_U = w c^T
    return S, blocks, labs


def comparison(S, blocks, scale=None):
    names = list(blocks)
    Bm = {n: np.array(blocks[n]).T for n in names}
    if scale is not None:
        for n in names:
            Bm[n] = Bm[n] * scale.get(n, 1.0)
    nb = len(names)
    C = np.zeros((nb, nb))
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            blk = Bm[a].T @ S @ Bm[b]
            if i == j:
                C[i, i] = np.linalg.eigvalsh((blk + blk.T) / 2)[0]
            else:
                C[i, j] = -np.linalg.norm(blk, 2)
    # Full matrix in these coordinates
    Ball = np.hstack([Bm[n] for n in names])
    Sf = Ball.T @ S @ Ball
    return names, C, np.linalg.eigvalsh((Sf + Sf.T) / 2)[0]


if __name__ == "__main__":
    rng = np.random.default_rng(7)
    for (k, m) in [(3, 1), (4, 1), (4, 2), (5, 1), (5, 2), (5, 3), (6, 2), (6, 3)]:
        worst = 9; worstw = None
        for trial in range(40):
            w = np.exp(rng.uniform(-2, 2, k))
            if trial % 4 == 1: w[:m] *= np.exp(rng.uniform(-5, -1))
            if trial % 4 == 2: w[m:] *= np.exp(rng.uniform(-5, -1))
            T = list(range(m))
            S, blocks, labs = reduced(w, T)
            names, C, lmin = comparison(S, blocks)
            cmin = np.linalg.eigvalsh(C)[0] / max(C.diagonal())
            # with blocks normalised so that diag = 1 (comparison matrix test is scale invariant per block)
            D = np.sqrt(np.abs(C.diagonal()))
            Cn = C / np.outer(D, D)
            cm = np.linalg.eigvalsh(Cn)[0]
            if cm < worst:
                worst = cm; worstw = (w ** 2 / (w ** 2).sum(), C.diagonal())
        print(f"k={k} m={m}: worst lambda_min of normalised comparison matrix {worst:.3f}; s/|s| = {np.round(worstw[0], 4)}, diag {np.round(worstw[1], 5)}", flush=True)
