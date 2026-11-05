"""Independent checks of the blocks used in the certificates (k = 4).

(1) Dense floating-point construction from the definitions: M = (I - 1/2 |Omega><Omega|)^{(x)4} on
    (C^3_row (x) C^3_col)^{(x)4}, V_p = p (x) I, D = V^* M V, R = M V D^{-1} V^* M, Y[(a,y'),(b,y)] = (M - R)[(b,y'),(a,y)],
    rescaled by mu_X = prod_t w_t^{N_t(X)}; compared with the exact ylib block at the same weights (3 weight vectors).
(2) The symbolic blocks of blocks_sym.pkl (s0 = 1), evaluated at rational points, equal the exact ylib blocks built
    over Q at those points (6 points, exact), including a point with s0 != 1 after the scaling s -> s/s0
    (entries homogeneous of degree n1 = j + 1 - |T|).
PASS/FAIL."""
import itertools, pickle, time, random
import numpy as np
import flint
from k4lib import block_states, pt_builder, BLOCKS
from ylib import Nv

t0 = time.time()
ok = True
k = 4
rng = np.random.default_rng(3)
Om = np.zeros(9)
Om[[0, 4, 8]] = 1
M1 = np.eye(9) - 0.5 * np.outer(Om, Om)
Mfull = M1
for _ in range(k - 1):
    Mfull = np.kron(Mfull, M1)


def gidx(a, y):
    i = 0
    for s in range(k):
        i = i * 9 + a[s] * 3 + y[s]
    return i


cols = list(itertools.product(range(3), repeat=k))
for trial in range(3):
    svals = [flint.fmpq(int(x), 13) for x in rng.integers(1, 40, size=k)]
    w = np.sqrt([float(x.p) / float(x.q) for x in svals])
    V = np.zeros((9 ** k, 3 ** k))
    for ci, z in enumerate(cols):
        for t in range(k):
            a = tuple(1 if r == t else 0 for r in range(k))
            V[gidx(a, z), ci] += w[t]
    MV = Mfull @ V
    Dm = V.T @ MV
    Dinv = np.linalg.inv(Dm)
    Bd = pt_builder(svals)
    for (j, T) in BLOCKS:
        S = block_states(j, T)
        idx = [gidx(b, yp) for (b, yp) in S]
        # Y[(a,y'),(b,y)] = W[(b,y'),(a,y)], W = M - R
        Yd = np.zeros((len(S), len(S)))
        for i1, (a, yp) in enumerate(S):
            for i2, (b, y) in enumerate(S):
                r1, r2 = gidx(b, yp), gidx(a, y)
                Yd[i1, i2] = Mfull[r1, r2] - MV[r1] @ Dinv @ MV[r2]
        mu = np.array([np.prod(w ** np.array(Nv(a, y))) for (a, y) in S])
        Yd = mu[:, None] * Yd * mu[None, :]
        Ye = np.array([[float(x.p) / float(x.q) for x in row] for row in Bd.matrix(S)])
        err = np.abs(Yd - Ye).max()
        ok &= err < 1e-9
        print(f'  dense vs exact, s = {[str(x) for x in svals]}, j = {j}, T = {T}: max deviation {err:.1e}', flush=True)
print(f'(1) done [{time.time() - t0:.1f}s]', flush=True)

D = pickle.load(open('blocks_sym.pkl', 'rb'))
ctx = flint.fmpq_mpoly_ctx.get(tuple(D['names']), 'lex')
random.seed(4)
for trial in range(6):
    s0 = flint.fmpq(1) if trial < 4 else flint.fmpq(random.randint(2, 9), random.randint(2, 9))
    svals = [s0] + [flint.fmpq(random.randint(1, 50), random.randint(1, 50)) for _ in range(3)]
    Bd = pt_builder(svals)
    sc = [x / s0 for x in svals[1:]]
    for (j, T) in BLOCKS:
        blk = D['blocks'][(j, T)]
        S = blk['states']
        Me = Bd.matrix(S)
        n1 = j + 1 - len(T)
        good = True
        for a in range(len(S)):
            for b in range(len(S)):
                nn, dd = blk['M'][a][b]
                v = ctx.from_dict(nn)(*sc) / ctx.from_dict(dd)(*sc)
                good &= v * s0 ** n1 == Me[a][b]
        ok &= good
        print(f'  symbolic (s0 = 1) vs exact point builder, s = {[str(x) for x in svals]}, j = {j}, T = {T}: {good}', flush=True)
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
