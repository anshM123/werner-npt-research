"""Check of Theorem 4 (tangent-class/README.md): if X has rank <= 2 and both range(X) and range(X^*) lie in the tangent
space T_e (e = |0...0>, plus all single excitations), then q_n(X) >= 2^{1-n} sum_{s != t} ||Y_st||^2, where Y_st is the
2 x 2 block of X between the excitations at copies s (rows) and t (columns).
q_n is computed directly from partial traces: q_n(X) = sum_J (-1/2)^{|J|} ||Tr_J X||^2.
(a) random rank-two X in the class, n = 2..5; (b) adversarial minimisation of (q - bound)/||X||^2 at n = 3, 4;
(c) the bound is attained: the ratio (q - bound)/||X||^2 gets close to 0."""
import numpy as np, itertools
from scipy.optimize import minimize
rng = np.random.default_rng(21)

def ptrace_norms(X, n):
    """q_n(X) by direct partial traces."""
    T = X.reshape([3] * (2 * n)); tot = 0.0
    for k in range(n + 1):
        for J in itertools.combinations(range(n), k):
            Y = T
            for j in sorted(J, reverse=True):
                Y = np.trace(Y, axis1=j, axis2=Y.ndim // 2 + j)
            tot += (-0.5) ** k * np.sum(np.abs(Y) ** 2)
    return tot

def embed(n):
    """indices of T_e inside (C^3)^n: e first, then (s, z) for s = 0..n-1, z = 1, 2."""
    idx = [0]
    for s in range(n):
        for z in (1, 2):
            digits = [0] * n; digits[s] = z; idx.append(int(np.ravel_multi_index(digits, [3] * n)))
    return idx

def X_from(p, n):
    d = 2 * n + 1; P = p[:4 * d].reshape(2, 2 * d); Q = p[4 * d:].reshape(2, 2 * d)
    Ur = P[0, :d] + 1j * P[0, d:]; Ur2 = P[1, :d] + 1j * P[1, d:]; Vr = Q[0, :d] + 1j * Q[0, d:]; Vr2 = Q[1, :d] + 1j * Q[1, d:]
    small = np.outer(Ur, Vr.conj()) + np.outer(Ur2, Vr2.conj())          # (2n+1) x (2n+1), rank <= 2
    idx = embed(n); X = np.zeros((3 ** n, 3 ** n), complex); X[np.ix_(idx, idx)] = small
    return X, small

def bound(small, n):
    tot = 0.0
    for s in range(n):
        for t in range(n):
            if s != t:
                tot += np.sum(np.abs(small[1 + 2 * s:3 + 2 * s, 1 + 2 * t:3 + 2 * t]) ** 2)
    return 2.0 ** (1 - n) * tot

ok = True; worst = np.inf
for n in (2, 3, 4, 5):
    for tr in range(8):
        X, small = X_from(rng.normal(size=8 * (2 * n + 1)), n)
        r = (ptrace_norms(X, n) - bound(small, n)) / np.sum(np.abs(X) ** 2); worst = min(worst, r)
print(f"(a) random rank-two X with range and co-range in T_e, n = 2..5: smallest (q - bound)/||X||^2 = {worst:+.3e}")
ok &= worst > -1e-12
adv = np.inf
for n in (3, 4):
    for st in range(4):
        f = lambda p: (lambda X_s: (ptrace_norms(X_s[0], n) - bound(X_s[1], n)) / np.sum(np.abs(X_s[0]) ** 2))(X_from(p, n))
        res = minimize(f, rng.normal(size=8 * (2 * n + 1)), method='L-BFGS-B', options=dict(maxiter=300))
        adv = min(adv, res.fun)
print(f"(b) adversarial minimum of (q - bound)/||X||^2 at n = 3, 4: {adv:+.3e}")
ok &= adv > -1e-9
print(f"(c) the bound is tight to optimiser accuracy: {adv < 1e-4}")
ok &= adv < 1e-3
assert ok
print("PASS")
