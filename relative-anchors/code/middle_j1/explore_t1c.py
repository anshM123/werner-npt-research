"""Exploration: dependence of S_P[r, r'] (|T| = 1 block, Schur complement onto P in rescaled coordinates, scaled so that
|s| = 1) on the configuration of the other copies."""
import numpy as np
from j1lib import J1, block_states, mu


def SP(x, t=0):
    w = np.sqrt(np.asarray(x, dtype=float))
    k = len(w)
    J = J1(list(w), F=float)
    st = block_states(k, (t,))
    M = np.array(J.matrix(st), dtype=float)
    m = np.array([mu(z, list(w)) for z in st])
    Mr = m[:, None] * M * m[None, :]
    P = [i for i, (a, y) in enumerate(st) if a is not None and a[0] == 1]
    E = [i for i in range(len(st)) if i not in P]
    A = Mr[np.ix_(P, P)]; B = Mr[np.ix_(P, E)]; C = Mr[np.ix_(E, E)]
    return A - B @ np.linalg.solve(C, B.T)


# fix sigma = x_t, x_1, x_2; vary the rest
sig, x1, x2 = 0.2, 0.15, 0.25
rest = 1 - sig - x1 - x2
for conf in ([rest], [rest / 2, rest / 2], [rest * 0.9, rest * 0.1], [rest / 4] * 4, [rest / 10] * 10):
    x = [sig, x1, x2] + list(conf)
    S = SP(x)
    print(len(x), 'S[1,2] =', S[1, 2], ' S[0,1] =', S[0, 1], ' S[0,2] =', S[0, 2], ' S[1,1] =', S[1, 1])
