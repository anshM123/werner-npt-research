"""Lemma R on larger and degenerate block operators: zero blocks, rank-deficient, real, nilpotent, permutation-like,
unequal block sizes are not allowed (blocks share K), d up to 6, m up to 8."""
import numpy as np
rng = np.random.default_rng(7)
tn = lambda M: np.linalg.svd(M, compute_uv=False).sum()
fn2 = lambda M: float(np.real(np.vdot(M, M)))


def Rgap(Z, d, m):
    bl = Z.reshape(d, m, d, m)
    lhs = sum(tn(bl[a, :, b, :]) ** 2 for a in range(d) for b in range(d)) + fn2(sum(bl[a, :, a, :] for a in range(d)))
    rhs = sum(fn2(bl[a, :, b, :]) for a in range(d) for b in range(d)) + tn(Z) ** 2
    return (lhs - rhs) / max(tn(Z) ** 2, 1e-300)


worst = -np.inf
for t in range(4000):
    d = int(rng.integers(2, 7)); m = int(rng.integers(1, 9)); D = d * m
    kind = t % 8
    Z = rng.normal(size=(D, D)) + 1j * rng.normal(size=(D, D))
    if kind == 1:   # zero out random blocks
        mask = rng.random((d, d)) < 0.5
        Z = (Z.reshape(d, m, d, m) * mask[:, None, :, None]).reshape(D, D)
    elif kind == 2:  # rank 1
        Z = np.outer(rng.normal(size=D) + 1j * rng.normal(size=D), rng.normal(size=D) + 1j * rng.normal(size=D))
    elif kind == 3:  # real
        Z = rng.normal(size=(D, D))
    elif kind == 4:  # nilpotent (strictly upper triangular)
        Z = np.triu(Z, 1)
    elif kind == 5:  # block permutation with unitary blocks
        P = np.eye(d)[rng.permutation(d)]
        Z = np.zeros((D, D), complex)
        for a in range(d):
            b = int(np.argmax(P[a]))
            Q = np.linalg.qr(rng.normal(size=(m, m)) + 1j * rng.normal(size=(m, m)))[0]
            Z[a * m:(a + 1) * m, b * m:(b + 1) * m] = Q * rng.uniform(0.1, 2)
    elif kind == 6:  # equality family + tiny noise
        v = rng.normal(size=d) + 1j * rng.normal(size=d); v /= np.linalg.norm(v)
        W = rng.normal(size=(m, m)) + 1j * rng.normal(size=(m, m))
        Z = np.kron(np.outer(v, v.conj()), W) + 1e-6 * Z
    elif kind == 7:  # Hermitian indefinite
        Z = Z + Z.conj().T
    worst = max(worst, Rgap(Z, d, m))
print('Lemma R, 4000 degenerate/large cases: max normalised violation =', worst)
