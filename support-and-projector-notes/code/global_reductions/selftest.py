import numpy as np
from qtools import q_dense, q_sparse_np, dense_to_sparse, rand_c
rng = np.random.default_rng(1)
for n in (1, 2, 3):
    X = rand_c((3 ** n, 3 ** n), rng)
    a = q_dense(X, n); b = q_sparse_np(*dense_to_sparse(X, n), )
    assert abs(a - b) < 1e-9 * max(1, abs(a)), (n, a, b)
    # subset formula directly for n = 2
    vs = [rand_c(3, rng) for _ in range(n)]
    v = vs[0]
    for u in vs[1:]:
        v = np.kron(v, u)
    P = np.outer(v, v.conj()); P /= np.trace(P).real
    # product rank one: q = prod (||p||^2 - 1/2 (tr p)^2) = prod (1 - 1/2) = 2^-n
    assert abs(q_dense(P, n) - 2.0 ** (-n)) < 1e-12
print('PASS')
