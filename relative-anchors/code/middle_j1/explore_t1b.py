"""Exploration: |T| = 1 block, Schur complement onto the P-states (1@r, 2@t) after eliminating R and Q, in rescaled
coordinates.  Is it a Laplacian (zero row sums, nonpositive off-diagonal)?  Other orders too."""
import numpy as np
from j1lib import J1, block_states, mu

rng = np.random.default_rng(4)


def block(w, t=0):
    k = len(w)
    J = J1(list(w), F=float)
    st = block_states(k, (t,))
    M = np.array(J.matrix(st), dtype=float)
    m = np.array([mu(z, list(w)) for z in st])
    return st, M, (m[:, None] * M * m[None, :])


def idx(st, kind):
    out = []
    for i, (a, y) in enumerate(st):
        if kind == 'R' and a is not None and a[0] == 2: out.append(i)
        if kind == 'P' and a is not None and a[0] == 1: out.append(i)
        if kind == 'Q' and a is None: out.append(i)
    return out


def schur(Mr, keep, elim):
    A = Mr[np.ix_(keep, keep)]; B = Mr[np.ix_(keep, elim)]; C = Mr[np.ix_(elim, elim)]
    return A - B @ np.linalg.solve(C, B.T), np.linalg.eigvalsh(C)[0]


res = {'P|RQ offmax': [], 'P|RQ rowsum': [], 'elim RQ mineig': [], 'P|Q offmax': [], 'P|R offmax': []}
for trial in range(400):
    k = int(rng.integers(3, 10))
    w = np.sqrt(rng.dirichlet(np.ones(k) * rng.choice([0.2, 1, 4])) + 1e-12)
    t = int(rng.integers(0, k))
    st, M, Mr = block(w, t)
    R, P, Q = idx(st, 'R'), idx(st, 'P'), idx(st, 'Q')
    S, me = schur(Mr, P, R + Q)
    off = S[~np.eye(len(P), dtype=bool)]
    res['P|RQ offmax'].append(off.max() / np.abs(np.diag(S)).max())
    res['P|RQ rowsum'].append(np.abs(S.sum(axis=1)).max())
    res['elim RQ mineig'].append(me)
    S2, _ = schur(Mr, P, Q)
    res['P|Q offmax'].append(S2[~np.eye(len(P), dtype=bool)].max())
    S3, _ = schur(Mr, P, R)
    res['P|R offmax'].append(S3[~np.eye(len(P), dtype=bool)].max())
for key, v in res.items():
    print(f'{key}: min {np.min(v):.4g}, max {np.max(v):.4g}')
