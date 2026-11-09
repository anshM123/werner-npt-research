"""Stabilizer domination on path graphs (1D cluster states) by transfer matrices.
Path with edge weights c_1..c_{N-1} in {1,2}: (Gamma a)_i = c_{i-1} a_{i-1} + c_i a_{i+1}.
S(d, v) = sum_a prod_i omega^{v_i a_i} h(a_i, (Gamma a)_i + d_i),  h(0,0) = 1 - 3t, else 1.
State after site i: vector over (a_i, a_{i+1}).  usage: python path_transfer.py mode N trials seed t
mode = check (compare with brute force for small N) | random (random c, d, v) | anneal (maximise the ratio)"""
import sys, itertools
import numpy as np

mode, N, trials, seed, t = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
rng = np.random.default_rng(seed)
lam = 1 - 3 * t; tau = t / (1 - t)
om = np.exp(2j * np.pi / 3)


def S_transfer(c, d, v):
    n = len(d)
    cc = [0] + list(c) + [0]              # cc[i] = weight of edge (i, i+1) in 1-based, cc[0] = cc[n] = 0
    # vec[(a_{i-1}, a_i)] before processing site i; start a_0 = 0, a_1 free
    vec = np.zeros((3, 3), complex)
    vec[0, :] = 1
    for i in range(1, n + 1):
        new = np.zeros((3, 3), complex)
        for p, q, r in itertools.product(range(3), repeat=3):
            if i == n and r != 0:
                continue                  # a_{n+1} = 0
            field = (cc[i - 1] * p + cc[i] * r + d[i - 1]) % 3
            hval = lam if (q == 0 and field == 0) else 1.0
            new[q, r] += vec[p, q] * om ** ((v[i - 1] * q) % 3) * hval
        vec = new
    return vec[:, 0].sum()


def S_brute(c, d, v):
    n = len(d)
    G = np.zeros((n, n), dtype=np.int64)
    for i in range(n - 1):
        G[i, i + 1] = G[i + 1, i] = c[i]
    A = np.array(list(itertools.product(range(3), repeat=n)), dtype=np.int64)
    w = ((A != 0) | ((A @ G.T + d) % 3 != 0)).sum(1)
    return np.sum(om ** ((A @ v) % 3) * lam ** (n - w))


def ratio(c, d, v):
    s0 = S_transfer(c, d, np.zeros_like(v)).real
    return abs(S_transfer(c, d, v)) / s0, s0


if mode == 'check':
    worst = 0
    for _ in range(trials):
        c = rng.integers(1, 3, size=N - 1); d = rng.integers(0, 3, size=N); v = rng.integers(0, 3, size=N)
        worst = max(worst, abs(S_transfer(c, d, v) - S_brute(c, d, v)))
    print(f"N={N}: transfer vs brute force, max |difference| = {worst:.2e}")
elif mode == 'random':
    best = 0; mins0 = np.inf
    for _ in range(trials):
        c = rng.integers(1, 3, size=N - 1); d = rng.integers(0, 3, size=N); v = rng.integers(0, 3, size=N)
        if not v.any():
            continue
        r, s0 = ratio(c, d, v); best = max(best, r); mins0 = min(mins0, s0)
    print(f"N={N} t={t}: random (c,d,v): max |S(d,v)|/S(d,0) = {best:.10f} (tau = {tau:.10f}); min S(d,0) = {mins0:.4e}")
else:
    best = 0
    for r_ in range(trials):
        c = rng.integers(1, 3, size=N - 1); d = rng.integers(0, 3, size=N); v = rng.integers(0, 3, size=N)
        if not v.any():
            v[0] = 1
        cur, _ = ratio(c, d, v)
        for it in range(400):
            c2, d2, v2 = c.copy(), d.copy(), v.copy()
            k = rng.integers(0, 3)
            if k == 0 and N > 1:
                c2[rng.integers(0, N - 1)] = rng.integers(1, 3)
            elif k == 1:
                d2[rng.integers(0, N)] = rng.integers(0, 3)
            else:
                v2[rng.integers(0, N)] = rng.integers(0, 3)
            if not v2.any():
                continue
            r2, _ = ratio(c2, d2, v2)
            if r2 >= cur:
                c, d, v, cur = c2, d2, v2, r2
        best = max(best, cur)
    print(f"N={N} t={t}: annealed max |S(d,v)|/S(d,0) = {best:.10f} (tau = {tau:.10f})")
