"""Exploration (floating point): the letter-1 block of (R_inf)_1, size (3k^2 - k)/2, rescaled coordinates.
States P_{r,c} = (1@r, 1@c), Q_{cc'} = (vac, 1@{c,c'}).  Kernel: p (x) |1@c> = indicator of the column group c.
(1) gap off the kernel (gauge-fixed smallest eigenvalue, and smallest nonzero eigenvalue of the block) against
    min_h s_h/|s|;
(2) Schur complement onto one column group {(r, c): r} after eliminating all other states: Laplacian?
(3) Schur complement onto the P-states after eliminating the vacuum states: sign pattern by type;
(the Schur complement onto a single row group is not defined: the rest is singular)."""
import numpy as np
from j1lib import J1, block_states, mu

rng = np.random.default_rng(12)


def block(x):
    w = np.sqrt(np.asarray(x, dtype=float))
    k = len(w)
    J = J1(list(w), F=float)
    st = block_states(k, ())
    M = np.array(J.matrix(st), dtype=float)
    m = np.array([mu(z, list(w)) for z in st])
    return st, m[:, None] * M * m[None, :], M


def schur(M, keep, elim):
    return M[np.ix_(keep, keep)] - M[np.ix_(keep, elim)] @ np.linalg.solve(M[np.ix_(elim, elim)], M[np.ix_(elim, keep)])


stats = {}
for trial in range(150):
    k = int(rng.integers(4, 8))   # k = 3 is (R_inf)_{k-2}: extra kernel vector (Proposition K2)
    x = rng.dirichlet(np.ones(k) * rng.choice([0.2, 1, 4])) + 1e-9
    x /= x.sum()
    st, Mr, M = block(x)
    P = [i for i, (a, y) in enumerate(st) if a is not None]
    Q = [i for i, (a, y) in enumerate(st) if a is None]
    col = lambda i: st[i][1].index(1)
    row = lambda i: st[i][0][1]
    # (1) gap: unrescaled block, smallest eigenvalue on the orthogonal complement of the kernel
    ev = np.linalg.eigvalsh(M)
    nz = ev[k:]   # k kernel vectors
    stats.setdefault('gap/min_x', []).append(nz[0] / x.min())
    stats.setdefault('kernel_eigs_max', []).append(np.abs(ev[:k]).max())
    # (2) column group c
    for c in range(k):
        keep = [i for i in P if col(i) == c]
        r0 = 0 if c != 0 else 1
        elim = [i for i in range(len(st)) if i not in keep and not (i in P and row(i) == r0)]   # gauge: drop (1@r0, 1@c')
        S = schur(Mr, keep, elim)
        off = S[~np.eye(k, dtype=bool)]
        stats.setdefault('colgroup offmax/diag', []).append(off.max() / np.diag(S).max())
        stats.setdefault('colgroup rowsum', []).append(np.abs(S.sum(1)).max() / np.diag(S).max())
    # (3) eliminate the vacuum only
    S = schur(Mr, P, Q)
    for a, i in enumerate(P):
        for b, j in enumerate(P):
            if a >= b:
                continue
            ri, ci, rj, cj = row(i), col(i), row(j), col(j)
            if ci == cj: typ = 'same col'
            elif ri == rj: typ = 'same row'
            elif ri == cj and ci == rj: typ = 'transposed'
            elif ri == ci or rj == cj: typ = 'one diagonal'
            else: typ = 'generic'
            stats.setdefault('P|Q ' + typ, []).append(S[a, b] / np.sqrt(S[a, a] * S[b, b]))
for key, v in stats.items():
    print(f'{key}: min {np.min(v):.4g}, max {np.max(v):.4g}')
