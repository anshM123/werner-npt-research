"""Stabilizer domination on CONNECTED graphs: max over connected Gamma, all shifts d and characters v != 0 of
|S(d,v)|/S(d,0), S(d,v) = sum_a omega^{v.a} (1-3t)^{N - w_d(a)}, w_d(a) = |supp a U supp(Gamma a + d)|.
Vectorised over d.  usage: python connected_domination.py N t mode count seed   (mode = all | sample)"""
import sys, itertools
import numpy as np
N, t, mode, count, seed = int(sys.argv[1]), float(sys.argv[2]), sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
lam = 1 - 3 * t; tau = t / (1 - t)
rng = np.random.default_rng(seed)
A = np.array(list(itertools.product(range(3), repeat=N)), dtype=np.int64)     # 3^N x N, also the list of d
pairs = [(i, j) for i in range(N) for j in range(i + 1, N)]
pw = lam ** (N - np.arange(N + 1))
nzA = (A != 0)
def connected(G):
    seen = {0}; stack = [0]
    while stack:
        i = stack.pop()
        for j in range(N):
            if G[i, j] and j not in seen:
                seen.add(j); stack.append(j)
    return len(seen) == N
def graphs():
    if mode == 'all':
        for e in itertools.product(range(3), repeat=len(pairs)):
            G = np.zeros((N, N), dtype=np.int64)
            for (i, j), x in zip(pairs, e):
                G[i, j] = G[j, i] = x
            yield G
    else:
        for _ in range(count):
            G = np.zeros((N, N), dtype=np.int64)
            for (i, j) in pairs:
                G[i, j] = G[j, i] = rng.integers(0, 3)
            yield G
worst = 0.0; arg = None; nconn = 0; minS0 = np.inf
for G in graphs():
    if not connected(G):
        continue
    nconn += 1
    GA = (A @ G.T) % 3                                        # Gamma a, rows a
    # field[d, a, i] = (Gamma a + d)_i
    fld = (GA[None, :, :] + A[:, None, :]) % 3
    w = (nzA[None, :, :] | (fld != 0)).sum(2)                 # (3^N d) x (3^N a)
    f = pw[w].reshape((3 ** N,) + (3,) * N)
    S = np.fft.fftn(f, axes=tuple(range(1, N + 1))).reshape(3 ** N, 3 ** N)
    S0 = S[:, 0].real
    minS0 = min(minS0, S0.min())
    r = np.abs(S[:, 1:]).max(1) / S0
    k = int(np.argmax(r))
    if r[k] > worst + 1e-12:
        worst = r[k]; arg = (G.copy(), A[k].copy())
print(f"N={N} t={t} ({mode}): {nconn} connected graphs, max |S(d,v)|/S(d,0) = {worst:.10f}, tau = {tau:.6f}, "
      f"ratio/tau = {worst / tau:.6f}, min S(d,0) = {minS0:.4e}")
print("  extremal graph:", arg[0].tolist(), " d:", arg[1].tolist())
