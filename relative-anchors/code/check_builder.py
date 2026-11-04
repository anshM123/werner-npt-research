"""Cross-check ylib.Builder (rescaled, exact) against an independent dense floating-point construction:
M = (I - 1/2 |Omega><Omega|)^{(x)k} on (C^3_row (x) C^3_col)^{(x)k}, V_p = p (x) I, D = V^* M V, R = M V D^{-1} V^* M,
Y[(a,y'),(b,y)] = (M - R)[(b,y'),(a,y)], then rescaled by mu_X = prod_t w_t^{N_t(X)}.  k = 2, 3, every j. PASS/FAIL."""
import itertools
import numpy as np
import flint
from ylib import Builder, point_ring, rinf_states, Nv

ok = True
rng = np.random.default_rng(1)
for k in (2, 3):
    svals = [flint.fmpq(int(x), 7) for x in rng.integers(1, 20, size=k)]
    w = np.sqrt([float(x.p) / float(x.q) for x in svals])
    # dense one-copy M on (row, col) with index row*3 + col
    Om = np.zeros(9); Om[[0, 4, 8]] = 1
    M1 = np.eye(9) - 0.5 * np.outer(Om, Om)
    M = M1
    for _ in range(k - 1):
        M = np.kron(M, M1)
    # global index: copies in order, each copy (row, col): idx = sum over copies
    def gidx(a, y):
        i = 0
        for s in range(k):
            i = i * 9 + a[s] * 3 + y[s]
        return i
    cols = list(itertools.product(range(3), repeat=k))
    p = {}
    for t in range(k):
        a = tuple(1 if r == t else 0 for r in range(k)); p[a] = w[t]
    V = np.zeros((9 ** k, 3 ** k))
    for ci, z in enumerate(cols):
        for a, wa in p.items():
            V[gidx(a, z), ci] += wa
    D = V.T @ M @ V
    R = M @ V @ np.linalg.solve(D, V.T @ M)
    W = M - R
    Bd = Builder(point_ring(svals))
    for j in range(-1, k + 1):
        S = rinf_states(k, j)
        Yd = np.array([[W[gidx(b, yp), gidx(a, y)] for (b, y) in S] for (a, yp) in S])
        mu = np.array([np.prod(w ** np.array(Nv(a, y))) for (a, y) in S])
        Yd = mu[:, None] * Yd * mu[None, :]
        Ye = np.array([[float(x.p) / float(x.q) for x in row] for row in Bd.matrix(S)]) if S else np.zeros((0, 0))
        err = np.abs(Yd - Ye).max() if S else 0.0
        ok &= err < 1e-10
        print(f'k={k} j={j}: {len(S)} states, max |dense - exact| = {err:.1e}, lmin = {np.linalg.eigvalsh(Ye).min() if S else 0:.2e}')
print('PASS' if ok else 'FAIL')
