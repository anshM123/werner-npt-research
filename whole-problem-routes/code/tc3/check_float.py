"""Floating-point check of the three-copy certificate from its definition (no orbit bookkeeping).

Decodes Z0, Z1 and L from cert/tc3_certificate.json.xz, builds the operators from diagram matrices, averages over
permutations of the copies, and compares at random real points
    lhs = f(x,y)|y|^2 W - t|x|^2|y|^4 W - t kap |x|^2|y|^2 sum_s <w_s|rho^y_s|w_s> W/|w_s|^2,   W = prod_s |w_s|^2,
    rhs = Avg_sigma [ <Psi|Z0 (x) I_w + Z1|Psi> + Tr(Lambda(x, y) G_0(y, w_0)) |w_1|^2 |w_2|^2 ] o sigma,
with f(x,y) = <y|P_{3/4}^{(x)3}(x x^T)|y> computed from the map itself.  Ends with PASS or FAIL."""
import itertools
import json
import lzma
import numpy as np
from exact_alg import build

CERT = json.load(lzma.open('cert/tc3_certificate.json.xz', 'rt', encoding='utf-8'))
t = 1 / 64
kap = 98 / 27
MAT = {k: np.array([M.astype(float) for M in build(k)['mats']]) for k in (2, 3, 4)}


def fl(sx):
    p, q = sx.split('/')
    return int(p) / int(q)


def decode(term):
    ks = CERT['terms'][term]
    algs = [CERT['algebras'][str(k)] for k in ks]
    cidx = [{tuple(c): i for i, c in enumerate(al['coords'])} for al in algs]
    z = np.zeros([len(al['coords']) for al in algs])
    for B in CERT['blocks']:
        if B['term'] != term:
            continue
        beta = B['beta']
        R = np.array(B['R'], dtype=float)
        Y = np.array([[fl(x) for x in row] for row in B['Y']])
        G = R @ Y @ R.T
        ms = [len(al['V'][str(j)][0]) for al, j in zip(algs, beta)]
        rows_ = list(itertools.product(*[range(m) for m in ms]))
        for Ri, rr in enumerate(rows_):
            for Ci, cc in enumerate(rows_):
                z[tuple(cidx[s][(beta[s], rr[s], cc[s])] for s in range(3))] += G[Ri, Ci]
    Psis = [np.array([[fl(x) for x in row] for row in al['Psi']]) for al in algs]
    return np.einsum('ai,bj,ck,ijk->abc', *Psis, z)


Z0d, Z1d, Ld = decode('Z0'), decode('Z1'), decode('L')
D2, D3, D4 = MAT[2], MAT[3], MAT[4]


def form(T, coeff, Dfirst):
    tot = 0.0
    for a0 in range(coeff.shape[0]):
        if not np.any(coeff[a0]):
            continue
        T0 = np.tensordot(Dfirst[a0], T, axes=([1], [0]))
        for a1 in range(15):
            if not np.any(coeff[a0, a1]):
                continue
            T01 = np.tensordot(D3[a1], T0, axes=([1], [1])).transpose(1, 0, 2)
            G = np.tensordot(T, T01, axes=([0, 1], [0, 1]))
            tot += np.einsum('akl,kl->a', D3, G) @ coeff[a0, a1]
    return tot


Lop = np.zeros((27 * 81, 27 * 81))
for a0 in range(15):
    for a1 in range(3):
        for a2 in range(3):
            if Ld[a0, a1, a2] != 0:
                Lop += Ld[a0, a1, a2] * np.kron(np.kron(D3[a0], D2[a1]), D2[a2])
Lr = Lop.reshape(3, 729, 3, 729)


def pmap(X):
    T = X.reshape([3] * 6)
    for s in range(3):
        T = np.moveaxis(T, [s, 3 + s], [0, 1])
        sh = T.shape
        Mi = T.reshape(3, 3, -1)
        tr = np.einsum('iir->r', Mi)
        T = np.moveaxis((-0.25 * (Mi + Mi.transpose(1, 0, 2)) + 0.75 * np.eye(3)[:, :, None] * tr).reshape(sh), [0, 1], [s, 3 + s])
    return T.reshape(27, 27)


def marg(v, s):
    V = np.moveaxis(v.reshape(3, 3, 3), s, 0).reshape(3, -1)
    return V @ V.T


def rhs(x, y, w):
    X, Y = x.reshape(3, 3, 3), y.reshape(3, 3, 3)
    tot = 0.0
    for perm in itertools.permutations(range(3)):
        Xp, Yp, wp = X.transpose(perm), Y.transpose(perm), [w[p] for p in perm]
        nw = [v @ v for v in wp]
        T3 = np.einsum('abc,def,ghi->adgbehcfi', Xp, Yp, Yp).reshape(27, 27, 27)
        tot += form(T3, Z0d, D3) * np.prod(nw)
        T4 = np.einsum('abc,def,ghi,j->adgjbehcfi', Xp, Yp, Yp, wp[0]).reshape(81, 27, 27)
        tot += form(T4, Z1d, D4) * nw[1] * nw[2]
        XY = np.einsum('abc,def->adbecf', Xp, Yp).reshape(-1)
        Lam = np.einsum('p,ipjq,q->ij', XY, Lr, XY)
        rho0 = marg(Yp.reshape(-1), 0)
        G0 = nw[0] * rho0 - (wp[0] @ rho0 @ wp[0]) * np.eye(3)
        tot += np.trace(Lam @ G0) * nw[1] * nw[2]
    return tot / 6


def lhs(x, y, w):
    nw = [v @ v for v in w]
    W = np.prod(nw)
    val = (y @ pmap(np.outer(x, x)) @ y) * (y @ y) * W - t * (x @ x) * (y @ y) ** 2 * W
    for s in range(3):
        val -= t * kap * (x @ x) * (y @ y) * (w[s] @ marg(y, s) @ w[s]) * W / nw[s]
    return val


rng = np.random.default_rng(29)
worst = 0.0
for trial in range(4):
    x, y = rng.standard_normal(27), rng.standard_normal(27)
    w = [rng.standard_normal(3) for _ in range(3)]
    a, b = lhs(x, y, w), rhs(x, y, w)
    worst = max(worst, abs(a - b) / max(1.0, abs(a)))
    print(f'random point {trial}: lhs {a:+.12e}  rhs {b:+.12e}', flush=True)
print(f'largest relative difference {worst:.1e}')
print('PASS' if worst < 1e-10 else 'FAIL')
