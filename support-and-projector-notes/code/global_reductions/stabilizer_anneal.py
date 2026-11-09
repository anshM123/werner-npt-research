"""Simulated annealing over qutrit stabilizer codes [[N,1]]_3 to minimise q_N(t, Pi)/||Pi||^2 over 2-dim code subspaces
Pi = P_code - |c><c| (c a logical stabilizer state).  Exact enumerator formula (validated against direct partial traces at N = 5):
    q = (4/3^N) [ W_S + W_l / 2 ],  W_X = sum_{E in X} (1-3t)^{N - wt E}.
Moves: random symplectic transvections u -> u + <u,v> v (v random, or supported on few sites).
usage: python stab_anneal.py N t steps restarts seed"""
import sys
import numpy as np

N, t, steps, restarts, seed = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
rng = np.random.default_rng(seed)
lam = 1 - 3 * t
powtab = lam ** (N - np.arange(N + 1))


def transvect(M, v):
    coef = (M[:, :N] @ v[N:] - M[:, N:] @ v[:N]) % 3
    return (M + np.outer(coef, v)) % 3


def span(gens):
    pts = np.zeros((1, 2 * N), dtype=np.int64)
    for g in gens:
        pts = np.concatenate([(pts + k * g) % 3 for k in range(3)])
    return pts


def weights(P):
    return ((P[:, :N] != 0) | (P[:, N:] != 0)).sum(1)


def value(M):
    S, LX, LZ = M[:N - 1], M[N - 1], M[N]
    Sp = span(S)
    WS = powtab[weights(Sp)].sum()
    best = np.inf
    for v in (LZ, LX, (LX + LZ) % 3, (LX + 2 * LZ) % 3):
        Wl = powtab[weights((Sp + v) % 3)].sum()
        best = min(best, 4 / 3 ** N * (WS + 0.5 * Wl) / 2)
    return best


def start():
    M = np.zeros((N + 1, 2 * N), dtype=np.int64)
    for i in range(N - 1):
        M[i, N + i] = 1
    M[N - 1, N - 1] = 1
    M[N, 2 * N - 1] = 1
    for _ in range(12 * N):
        M = transvect(M, rng.integers(0, 3, size=2 * N))
    return M


def move(M):
    v = np.zeros(2 * N, dtype=np.int64)
    k = rng.integers(1, 4)
    sites = rng.choice(N, size=k, replace=False)
    for s in sites:
        v[s] = rng.integers(0, 3); v[N + s] = rng.integers(0, 3)
    if not v.any():
        v[rng.integers(0, 2 * N)] = 1
    return transvect(M, v)


overall = (np.inf, None)
for r in range(restarts):
    M = start(); cur = value(M); best = (cur, M.copy())
    T0 = 0.05
    for it in range(steps):
        T = T0 * (1 - it / steps) + 1e-5
        M2 = move(M); v2 = value(M2)
        if v2 < cur or rng.random() < np.exp(-(v2 - cur) / T):
            M, cur = M2, v2
            if cur < best[0]:
                best = (cur, M.copy())
    print(f"N={N} t={t} restart {r}: best q/||Pi||^2 = {best[0]:.6f}", flush=True)
    if best[0] < overall[0]:
        overall = best
pass
print(f"N={N} t={t}: overall best q/||Pi||^2 over annealed code subspaces = {overall[0]:.6f}")
