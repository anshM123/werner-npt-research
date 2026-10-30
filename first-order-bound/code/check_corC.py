"""Check Corollary C: S_m(Z) <= C(n-1,m)||Z||_2^2 + C(n-1,m-1)||Z||_1^2 for random Z, all m."""
import itertools
from math import comb
import numpy as np
rng = np.random.default_rng(99)
worst = -np.inf
for trial in range(600):
    d = int(rng.integers(2, 4)); n = int(rng.integers(2, 6 if d == 2 else 4))
    D = d ** n
    if trial % 3 == 0:
        Z = rng.normal(size=(D, D)) + 1j * rng.normal(size=(D, D))
    elif trial % 3 == 1:
        r = int(rng.integers(1, 3)); Z = (rng.normal(size=(D, r)) + 1j * rng.normal(size=(D, r))) @ (rng.normal(size=(r, D)) + 1j * rng.normal(size=(r, D)))
    else:   # near equality: P (x) |c><c|
        P = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)); P = P @ P.conj().T
        c = np.ones(1)
        for _ in range(n - 1):
            u = rng.normal(size=d) + 1j * rng.normal(size=d); c = np.kron(c, u / np.linalg.norm(u))
        Z = np.kron(P, np.outer(c, c.conj())) + 1e-3 * (rng.normal(size=(D, D)) + 1j * rng.normal(size=(D, D)))
    T = Z.reshape([d] * (2 * n))
    n2 = np.real(np.vdot(Z, Z)); n1 = np.linalg.svd(Z, compute_uv=False).sum()
    for m in range(1, n + 1):
        S = 0.0
        for J in itertools.combinations(range(n), m):
            cur = T; k = n
            for i in sorted(J, reverse=True):
                cur = np.trace(cur, axis1=i, axis2=k + i); k -= 1
            S += np.real(np.vdot(cur, cur))
        bound = comb(n - 1, m) * n2 + comb(n - 1, m - 1) * n1 ** 2
        worst = max(worst, (S - bound) / n1 ** 2)
print('Corollary C: max (S_m - bound)/||Z||_1^2 =', worst)
