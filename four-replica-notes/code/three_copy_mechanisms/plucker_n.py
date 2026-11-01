"""n-copy Plucker inequality (PSD endpoint with margin) and cut-wise strengthenings.
For orthonormal u, v in (C^d)^{(x)n}:  p_S = ||Pi_S (u (x) v)||^2, a_S(u) = ||Pi_S (u (x) u)||^2 (Pi_S: antisymmetric on S).
  R_n = sum_{|S| odd} (3^|S| - 3) p_S - sum_{|S| even} (3^|S| - 1) p_S
  target:  R_n <= sqrt(Q(u) Q(v)),  Q(u) = sum_{|S| even} (3^|S| - 1) a_S(u)
  variant 'pattern':  R_n <= sum_{|S| even >= 2} (3^|S| - 1) sqrt(a_S(u) a_S(v))
  variant 'size':     R_n <= sum_{k even >= 2} (3^k - 1) sqrt(A_k(u) A_k(v)),  A_k = sum_{|S| = k} a_S
Maximise LHS - RHS over u, v (Gram-Schmidt inside).  Pattern weights via local swap expansions."""
import sys, itertools, numpy as np
from scipy.optimize import minimize
d, n, variant, starts, seed = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
D = d ** n
subsets = [S for m in range(n + 1) for S in itertools.combinations(range(n), m)]
def pattern_weights(x, y):
    # w_S = ||Pi_S (x (x) y)||^2 with Pi_S = prod_{i in S} (I - F_i)/2 prod_{i notin S} (I + F_i)/2
    T = np.outer(x, y).reshape((d,) * (2 * n))
    w = {}
    # compute <x y| F_A |x y> for all A, then Mobius-type transform
    F = {}
    for A in subsets:
        TA = T
        for k in A:
            TA = np.swapaxes(TA, k, n + k)
        F[A] = np.vdot(T, TA).real
    for S in subsets:
        tot = 0.0
        for A in subsets:
            sign = (-1) ** len(set(A) & set(S))
            tot += sign * F[A]
        w[S] = tot / 2 ** n
    return w
def obj(v):
    v = v.reshape(4, D); u = v[0] + 1j * v[1]; w2 = v[2] + 1j * v[3]
    u = u / np.linalg.norm(u); w2 = w2 - np.vdot(u, w2) * u; w2 = w2 / np.linalg.norm(w2)
    p = pattern_weights(u, w2); au = pattern_weights(u, u); av = pattern_weights(w2, w2)
    R = sum((3 ** len(S) - 3) * p[S] for S in subsets if len(S) % 2 == 1) - sum((3 ** len(S) - 1) * p[S] for S in subsets if len(S) % 2 == 0)
    if variant == 'target':
        Qu = sum((3 ** len(S) - 1) * au[S] for S in subsets if len(S) % 2 == 0)
        Qv = sum((3 ** len(S) - 1) * av[S] for S in subsets if len(S) % 2 == 0)
        rhs = np.sqrt(max(Qu, 0) * max(Qv, 0))
    elif variant == 'pattern':
        rhs = sum((3 ** len(S) - 1) * np.sqrt(max(au[S], 0) * max(av[S], 0)) for S in subsets if len(S) % 2 == 0 and len(S) >= 2)
    else:
        rhs = 0.0
        for k in range(2, n + 1, 2):
            Au = sum(au[S] for S in subsets if len(S) == k); Av = sum(av[S] for S in subsets if len(S) == k)
            rhs += (3 ** k - 1) * np.sqrt(max(Au, 0) * max(Av, 0))
    return -(R - rhs)
rng = np.random.default_rng(seed); best = -np.inf
for s in range(starts):
    res = minimize(obj, rng.standard_normal(4 * D), method='L-BFGS-B', options={'maxiter': 3000})
    best = max(best, -res.fun)
    print(f"start {s}: max(LHS-RHS) {-res.fun:+.5f}", flush=True)
print(f"d={d} n={n} variant={variant}: max LHS-RHS = {best:+.6f}")
