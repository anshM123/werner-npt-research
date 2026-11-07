"""Numerics for Conjecture PM (the least value is attained on planes).

For parameters a = (a_1, ..., a_n) with a_s >= 1/2, the planar candidate value is
    m(a) = min{ prod_s (a_s - 1/2),  min over |S| >= 2 of 2^-|S| [prod_{s in S} (2 a_s - 1/2) - 1/2] prod_{s not in S} a_s },
the value of aligned real products, and of the vectors of Theorem 7.1 on the copies in S with orthogonal real products at
the others.  For each parameter vector the script minimises <y| (x)_s P_{a_s}(x x^T) |y> over real unit x, y by
alternating minimisation from random starts and from 'thick' starts (one-copy marginals bounded away from rank 2), with
no planar seeds, and reports the least value found, m(a), their difference, and the least one-copy marginal eigenvalue
of the minimiser.  A value below m(a) by more than rounding would contradict the conjecture.

usage: python planar_minimum.py n starts seed
"""
import itertools, sys
import numpy as np

d = 3
n = int(sys.argv[1]); starts = int(sys.argv[2]); rng = np.random.default_rng(int(sys.argv[3]))


def apply_map(a_list, X):
    T = X.reshape([d] * (2 * n)).astype(float)
    for s, a in enumerate(a_list):
        T = np.moveaxis(T, [s, n + s], [0, 1]); sh = T.shape
        Mi = T.reshape(3, 3, -1); out = np.zeros_like(Mi); tr = np.einsum('iir->r', Mi)
        for k in range(3): out[k, k] += a * tr
        out -= 0.25 * (Mi + Mi.transpose(1, 0, 2))
        T = out.reshape(sh); T = np.moveaxis(T, [0, 1], [s, n + s])
    return T.reshape(d**n, d**n)


def planar_min(a_list):
    cands = [np.prod([a - 0.5 for a in a_list])]
    for r in range(2, n + 1):
        for S in itertools.combinations(range(n), r):
            cands.append(2.0**-r * (np.prod([2 * a_list[s] - 0.5 for s in S]) - 0.5)
                         * np.prod([a_list[s] for s in range(n) if s not in S]))
    return min(cands)


def thick_start():
    # a random state whose one-copy marginals all have least eigenvalue >= 0.05 (rejection sampling on mixtures)
    while True:
        x = rng.standard_normal(d**n); x /= np.linalg.norm(x)
        X = x.reshape([d] * n)
        if all(np.linalg.eigvalsh(np.moveaxis(X, s, 0).reshape(3, -1) @ np.moveaxis(X, s, 0).reshape(3, -1).T)[0] >= 0.05
               for s in range(n)):
            return x


if n == 3:
    tests = [[0.5] * 3, [0.55] * 3, [0.6] * 3, [0.64] * 3, [0.6469] * 3, [0.7] * 3, [0.75] * 3, [0.5, 0.6, 0.7],
             [0.5, 0.75, 0.75], [0.55, 0.65, 0.8], [0.6, 0.6, 1.0]]
else:
    tests = [[0.5] * 4, [0.6] * 4, [0.62] * 4, [0.66] * 4, [0.75] * 4, [0.52, 0.56, 0.6, 0.64], [0.5, 0.5, 0.5, 0.75],
             [0.55, 0.6, 0.7, 0.75]]
worst = -np.inf
for a_list in tests:
    best = np.inf; bestx = None
    for st in range(starts):
        x = thick_start() if st % 2 else rng.standard_normal(d**n)
        x /= np.linalg.norm(x)
        for it in range(80):
            W = apply_map(a_list, np.outer(x, x)); w, V = np.linalg.eigh((W + W.T) / 2); y = V[:, 0]
            W2 = apply_map(a_list, np.outer(y, y)); w2, V2 = np.linalg.eigh((W2 + W2.T) / 2); x = V2[:, 0]
        if w2[0] < best: best = w2[0]; bestx = x.copy()
    m = planar_min(a_list)
    X = bestx.reshape([d] * n)
    lam = [np.linalg.eigvalsh(np.moveaxis(X, s, 0).reshape(3, -1) @ np.moveaxis(X, s, 0).reshape(3, -1).T)[0]
           for s in range(n)]
    worst = max(worst, m - best)
    print(f"a = {a_list}: least value found {best:+.6e}; planar value m(a) = {m:+.6e}; found - m = {best - m:+.1e}; "
          f"least one-copy marginal eigenvalues of the minimiser {np.round(lam, 4)}", flush=True)
print(f"n = {n}: largest amount by which a found value lies below m(a): {worst:.1e} (rounding level means no violation)")
