"""W_r = (n+r-1) I - sum_i P+_{A_i B_i} (x) I on (C^d (x) C^d)^{(x)n}: >= 0 on Schmidt rank <= r (A|B), negative on
the rank-(r+1) vector vec(I_U (x) |c><c|), dim U = r+1. Vectors phi = vec(M) with M: B -> A of rank r."""
import numpy as np
rng = np.random.default_rng(1)


def expval(M, d, n, r):
    # <vec M| W_r |vec M> = (n+r-1)||M||^2 - sum_i ||Tr_i M||^2
    T = M.reshape([d] * (2 * n)); val = (n + r - 1) * np.real(np.vdot(M, M))
    for i in range(n):
        Ti = np.trace(T, axis1=i, axis2=n + i); val -= np.real(np.vdot(Ti, Ti))
    return val


worst = np.inf
for t in range(3000):
    d = int(rng.integers(2, 4)); n = int(rng.integers(1, 5 if d == 2 else 4)); r = int(rng.integers(1, d + 1))
    D = d ** n
    M = (rng.normal(size=(D, r)) + 1j * rng.normal(size=(D, r))) @ (rng.normal(size=(r, D)) + 1j * rng.normal(size=(r, D)))
    worst = min(worst, expval(M, d, n, r) / np.real(np.vdot(M, M)))
print('min <W_r> / ||phi||^2 over random Schmidt-rank-r vectors:', worst)
# sharpness: rank r+1 example
for d, n, r in ((3, 3, 1), (3, 3, 2), (4, 2, 3)):
    U = np.diag([1.0] * (r + 1) + [0.0] * (d - r - 1))
    c = np.zeros(d); c[0] = 1
    M = U
    for _ in range(n - 1):
        M = np.kron(M, np.outer(c, c))
    print(f'd={d} n={n} r={r}: <W_r> on rank-{r+1} example = {expval(M, d, n, r):+.3f}  (expected -(r+1) = {-(r+1)})')
