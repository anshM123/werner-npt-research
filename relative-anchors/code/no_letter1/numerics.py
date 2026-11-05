"""Floating-point illustration (no claim rests on it): the margins of the block and of the two inequalities over
random configurations, |T| = 2, ..., 14, |U| = 0, ..., 40.  Reported: min eigenvalue of M, max_i<>l F_il / (a_i a_l),
and min over pairs of (a_l F_ii + (x_T - a_i) F_il)/(a_l (c - a_i)), c = 2 - x_T (the corner scale)."""
import numpy as np, itertools
rng = np.random.default_rng(20261105)


def block(a, xU):
    n = len(a); xT = a.sum(); c = 2 - xT
    Ph = lambda d: np.sum(xU / (c - 2 * xU + d))
    M = 2 * np.eye(n + 1)
    for i, l in itertools.combinations(range(n), 2):
        Gi = a[i] / (c - a[i] + a[l]); Gl = a[l] / (c - a[l] + a[i])
        M[i, l] = M[l, i] = (1 + Gi) * (1 + Gl) / (2 + Gi + Gl + Ph(a[i] + a[l]))
    for i in range(n):
        g = a[i] / (c - a[i])
        M[i, n] = M[n, i] = -2 * (1 + g) / (2 + g + Ph(a[i]))
    M[n, n] = 2 / (1 + Ph(0) / 2)
    return M


worst = {}
for trial in range(6000):
    n = int(rng.integers(2, 15)); N = int(rng.choice([0, 1, 2, 3, 5, 10, 40]))
    kind = rng.integers(0, 4)
    if kind == 0:
        w = rng.random(n + N)
    elif kind == 1:
        w = np.exp(rng.normal(0, 3, n + N))
    elif kind == 2:
        w = np.ones(n + N); w[0] = 10 ** rng.uniform(1, 5)
    else:
        w = np.concatenate([np.exp(rng.normal(0, 1, n)), 10 ** rng.uniform(-6, 2) * rng.random(N)])
    w = w / w.sum(); a, xU = w[:n], w[n:]
    M = block(a, xU)
    ev = np.linalg.eigvalsh(M).min()
    F = M[:n, :n] - np.outer(M[:n, n], M[n, :n]) / M[n, n]
    xT = a.sum(); c = 2 - xT
    off = max(F[i, l] / (a[i] * a[l]) for i in range(n) for l in range(n) if i != l)
    row = min((a[l] * F[i, i] + (xT - a[i]) * F[i, l]) / (a[l] * (c - a[i])) for i in range(n) for l in range(n) if i != l)
    key = (n > 5, N > 0)
    wv = worst.setdefault(key, [1e9, -1e9, 1e9])
    wv[0] = min(wv[0], ev); wv[1] = max(wv[1], off); wv[2] = min(wv[2], row)
for key in sorted(worst):
    wv = worst[key]
    print(f'|T| {"> 5" if key[0] else "<= 5"}, U {"nonempty" if key[1] else "empty"}: min eig M {wv[0]:.3e}; '
          f'max F_il/(a_i a_l) {wv[1]:.4f}; min scaled pair-row margin {wv[2]:.4f}')
print('numerical illustration done')
