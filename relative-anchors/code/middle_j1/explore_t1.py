"""Exploration (floating point): the type I block T = {t} of (R_inf)_1, size 3k - 1.
States: R_c = (2@t, 1@c) (c in [k]), P_r = (1@r, 2@t) (r in [k]), Q_c = (vac, 2@t 1@c) (c != t).
Rescaled coordinates (ylib mu): R_c, Q_c by w_c, P_r by w_r.  Gauge: delete P_{r0}, r0 != t.
Look at: smallest eigenvalue off the kernel; Schur complements onto R, onto (R, Q), sign patterns."""
import numpy as np
from j1lib import J1, block_states, mu

rng = np.random.default_rng(3)


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


stats = {}
for trial in range(300):
    k = int(rng.integers(3, 9))
    w = np.sqrt(rng.dirichlet(np.ones(k) * rng.choice([0.3, 1, 3])) + 1e-9)
    t = 0
    st, M, Mr = block(w, t)
    R, P, Q = idx(st, 'R'), idx(st, 'P'), idx(st, 'Q')
    r0 = 1
    Pk = [i for i in P if st[i][0][1] != r0]
    keep = R + Pk + Q
    G = Mr[np.ix_(keep, keep)]
    ev = np.linalg.eigvalsh(G)
    stats.setdefault('mineig_gauge', []).append(ev[0])
    # Schur complement onto R after eliminating Pk and Q
    nR = len(R)
    Arr = G[:nR, :nR]; Arx = G[:nR, nR:]; Axx = G[nR:, nR:]
    exx = np.linalg.eigvalsh(Axx)[0]
    stats.setdefault('mineig_PQ', []).append(exx)
    S = Arr - Arx @ np.linalg.solve(Axx, Arx.T)
    off = S[~np.eye(nR, dtype=bool)]
    stats.setdefault('S_offdiag_max', []).append(off.max())
    stats.setdefault('S_offdiag_min', []).append(off.min())
    stats.setdefault('S_mineig', []).append(np.linalg.eigvalsh(S)[0])
    # Schur complement onto (R, Q) after eliminating Pk
    RQ = list(range(nR)) + list(range(nR + len(Pk), len(keep)))
    PP = list(range(nR, nR + len(Pk)))
    S2 = G[np.ix_(RQ, RQ)] - G[np.ix_(RQ, PP)] @ np.linalg.solve(G[np.ix_(PP, PP)], G[np.ix_(PP, RQ)])
    off2 = S2[~np.eye(len(RQ), dtype=bool)]
    stats.setdefault('S2_offdiag_max', []).append(off2.max())
    # sign pattern of the raw rescaled block
    Rb = Mr[np.ix_(R, R)]; offR = Rb[~np.eye(nR, dtype=bool)]
    stats.setdefault('RR_offdiag_max', []).append(offR.max())
    stats.setdefault('RR_offdiag_min', []).append(offR.min())
for key, v in stats.items():
    print(f'{key}: min {np.min(v):.4g}, max {np.max(v):.4g}')
