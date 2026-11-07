"""Numerical check for Corollary 7.10 (whole-problem-routes/README.md): at a = 3/4, how low can <y|P^n(xx^T)|y> go on inputs x whose one-copy marginals all have rank 3?
Projected gradient search: minimise the least output eigenvalue subject to lambda_min(rho_s) >= tau at every copy
(projected gradient steps, halved until the constraint holds; starts near the state sum_k |k, k+1, k+2, ...> whose
marginals are all I/3).  Reports the least value found for each tau.  usage: python thick_inputs.py n"""
import numpy as np, sys
n = int(sys.argv[1]); a = 0.75; rng = np.random.default_rng(5 + n)
d = 3; D = d ** n
e = np.eye(3); M = np.zeros((9, 9))
for i in range(3):
    for j in range(3):
        E = np.outer(e[i], e[j]); M[:, 3 * i + j] = (a * np.trace(E) * np.eye(3) - (E + E.T) / 4).reshape(9)
M4 = M.reshape(3, 3, 3, 3)
def Phi(X):
    T = X.reshape([3] * (2 * n))
    for c in range(n):
        T = np.moveaxis(np.tensordot(M4, T, axes=([2, 3], [c, n + c])), [0, 1], [c, n + c])
    return T.reshape(D, D)
def marg_min(x):
    T = x.reshape([3] * n); out = []
    for c in range(n):
        A = np.moveaxis(T, c, 0).reshape(3, -1); out.append(np.linalg.eigvalsh(A @ A.T)[0])
    return min(out)
def value(x):
    return np.linalg.eigvalsh(Phi(np.outer(x, x)))[0]
# projected search: gradient of lambda_min(Phi(xx^T)) w.r.t. x is 2 Phi*(yy^T) x with y the bottom eigenvector
for tau in (0.25, 0.15, 0.08, 0.04, 0.02):
    best = np.inf; bestm = None
    for s in range(12):
        # thick start: mix a random vector with a state whose marginals are all I/3 (a product of maximally entangled pairs
        # is not available for odd n, so use the 'shifted' state sum_k |k, k+1, k+2, ...> / sqrt3 plus noise)
        base = np.zeros(D)
        for k in range(3):
            base[int(''.join(str((k + c) % 3) for c in range(n)), 3)] = 1
        base /= np.linalg.norm(base)
        x = base + 0.15 * rng.normal(size=D) / np.sqrt(D); x /= np.linalg.norm(x)
        if marg_min(x) < tau:
            continue
        for it in range(400):
            w, V = np.linalg.eigh(Phi(np.outer(x, x))); y = V[:, 0]
            g = 2 * Phi(np.outer(y, y)) @ x           # P is self-adjoint
            xn = x - 0.2 * g; xn /= np.linalg.norm(xn)
            # keep the marginal constraint: if violated, mix back toward x
            t = 1.0
            while marg_min(xn) < tau and t > 1e-4:
                t /= 2; xn = x - 0.2 * t * g; xn /= np.linalg.norm(xn)
            if marg_min(xn) < tau:
                break
            x = xn
        v = value(x)
        if v < best:
            best, bestm = v, marg_min(x)
    if bestm is None:
        print(f'tau = {tau}: no admissible start'); continue
    print(f'n = {n}, a = 3/4, marginals >= {tau}: least value found {best:+.4e} (min marginal eigenvalue {bestm:.3f}); '
          f'product value (1/4)^n = {0.25 ** n:.4e}', flush=True)
