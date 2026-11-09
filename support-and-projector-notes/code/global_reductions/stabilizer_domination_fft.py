"""Stabilizer domination test via FFT over Z_3^N: for each shift d, f_d(a) = lam^{N - w_d(a)} on F_3^N and
S(d, .) = its Fourier transform (computed with numpy.fft.fftn on shape (3,)*N).  Reports max |S(d,v)|/S(d,0), v != 0.
usage: python fourier_dom_fft.py N t graphs seed"""
import sys, itertools
import numpy as np

N, t, graphs, seed = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
lam = 1 - 3 * t
rng = np.random.default_rng(seed)
A = np.array(list(itertools.product(range(3), repeat=N)), dtype=np.int64)
pw = lam ** (N - np.arange(N + 1))
worst = 0.0; ties = 0
for g in range(graphs):
    M = np.triu(rng.integers(0, 3, size=(N, N)), 1); G = (M + M.T) % 3
    GA = (A @ G.T) % 3
    nzA = (A != 0)
    for d in A:
        w = (nzA | ((GA + d) % 3 != 0)).sum(1)
        f = pw[w].reshape((3,) * N)
        # numpy fft uses exp(-2 pi i k n / 3); |S| is the same for v and -v, so the sign convention is immaterial
        S = np.fft.fftn(f).ravel()
        S0 = S[0].real
        r = np.abs(S[1:]).max() / S0
        worst = max(worst, r)
        ties += int((np.abs(np.abs(S[1:]) - S0) < 1e-9 * S0).sum())
    print(f"graph {g}: running max ratio {worst:.10f}", flush=True)
print(f"N={N} t={t}: {graphs} random graphs, all d and v: max |S(d,v)|/S(d,0) = {worst:.10f}; equality cases {ties}")
