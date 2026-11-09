import itertools, numpy as np
rng = np.random.default_rng(4)
om = np.exp(2j * np.pi / 3)
def S(G, d, v, t):
    N = len(d)
    if N == 0:
        return 1.0 + 0j
    A = np.array(list(itertools.product(range(3), repeat=N)), dtype=np.int64)
    w = ((A != 0) | ((A @ G.T + d) % 3 != 0)).sum(1)
    return np.sum(om ** ((A @ v) % 3) * (1 - 3 * t) ** (N - w))
worst = 0
for trial in range(60):
    N = int(rng.integers(1, 6)); t = rng.uniform(0.2, 0.5)
    M = np.triu(rng.integers(0, 3, size=(N, N)), 1); G = (M + M.T) % 3
    d = rng.integers(0, 3, size=N); v = rng.integers(0, 3, size=N)
    Gp, gam = G[:-1, :-1], G[-1, :-1]
    dp, dn, vp, vn = d[:-1], d[-1], v[:-1], v[-1]
    rec = sum(om ** (m * vn) * S(Gp, (dp + m * gam) % 3, vp, t) - t * om ** (m * dn) * S(Gp, dp, (vp + m * gam) % 3, t) for m in range(3))
    worst = max(worst, abs(rec - S(G, d, v, t)) / max(1, abs(S(G, d, v, t))))
print("vertex-deletion recursion: max relative error over 60 random (graph, d, v, t):", worst)
