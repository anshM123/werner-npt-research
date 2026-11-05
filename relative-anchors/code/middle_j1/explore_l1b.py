"""Exploration (floating point): sequential column-group elimination in the letter-1 block of (R_inf)_1 (k >= 4).
Rescaled coordinates.  Stage 0: eliminate the vacuum states.  Then, for the column groups c_1, c_2, ... in a random
order: the diagonal block of the current Schur complement on the group c_m should be a Laplacian (zero row sums,
nonpositive off-diagonal); then eliminate the group in the gauge (delete one coordinate (1@r0, 1@c_m)).
Also the same for the |T| = 1 block (one column group, plus the R and Q states) as a sanity check."""
import numpy as np
from j1lib import J1, block_states, mu

rng = np.random.default_rng(21)
worst_off = -1e9; worst_sum = 0; worst_piv = 1e9; n = 0
for trial in range(200):
    k = int(rng.integers(4, 8))
    x = rng.dirichlet(np.ones(k) * rng.choice([0.3, 1, 4])) + 1e-7
    x /= x.sum()
    w = np.sqrt(x)
    J = J1(list(w), F=float)
    st = block_states(k, ())
    M = np.array(J.matrix(st), dtype=float)
    m = np.array([mu(z, list(w)) for z in st])
    Mr = m[:, None] * M * m[None, :]
    P = [i for i, (a, y) in enumerate(st) if a is not None]
    Q = [i for i, (a, y) in enumerate(st) if a is None]
    col = {i: st[i][1].index(1) for i in P}
    row = {i: st[i][0][1] for i in P}
    # stage 0
    cur = P[:]
    S = Mr[np.ix_(P, P)] - Mr[np.ix_(P, Q)] @ np.linalg.solve(Mr[np.ix_(Q, Q)], Mr[np.ix_(Q, P)])
    order = list(rng.permutation(k))
    for c in order:
        g = [a for a, i in enumerate(cur) if col[i] == c]
        D = S[np.ix_(g, g)]
        off = D[~np.eye(len(g), dtype=bool)]
        sc = np.abs(np.diag(D)).max()
        worst_off = max(worst_off, off.max() / sc)
        worst_sum = max(worst_sum, np.abs(D.sum(1)).max() / sc)
        n += 1
        # eliminate the group in the gauge: keep all of it except one coordinate (row r0), which is dropped
        r0 = row[cur[g[0]]]
        el = g[1:]
        dropped = g[0]
        piv = np.linalg.eigvalsh(S[np.ix_(el, el)])[0]
        worst_piv = min(worst_piv, piv / sc)
        rest = [a for a in range(len(cur)) if a not in g]
        S = S[np.ix_(rest, rest)] - S[np.ix_(rest, el)] @ np.linalg.solve(S[np.ix_(el, el)], S[np.ix_(el, rest)])
        cur = [cur[a] for a in rest]
print(f'{n} stages: max off-diagonal/diag = {worst_off:.3g}, max |row sum|/diag = {worst_sum:.3g}, '
      f'min eigenvalue of the eliminated gauge-fixed group / diag = {worst_piv:.3g}')
