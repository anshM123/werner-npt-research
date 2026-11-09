"""Fourier domination for Pauli-shifted graph states: for symmetric Gamma over F_3, shift d and character v,
    S(d, v) = sum_a omega^{v.a} lam^{N - w_d(a)},  w_d(a) = |supp a  U  supp(Gamma a + d)|,
EP on the class G_N = {lam1 Z^x1 G Z^-y1 + lam2 Z^x2 G Z^-y2} is equivalent to |S(d, v)| <= S(d, 0) for all d, v.
Report max over (Gamma, d, v != 0) of |S(d,v)|/S(d,0) (exhaustive over Gamma for N <= 4, sampled otherwise).
usage: python fourier_dom.py N t mode count seed"""
import sys, itertools
import numpy as np

N, t, mode, count, seed = int(sys.argv[1]), float(sys.argv[2]), sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
lam = 1 - 3 * t
rng = np.random.default_rng(seed)
A = np.array(list(itertools.product(range(3), repeat=N)), dtype=np.int64)
pairs = [(i, j) for i in range(N) for j in range(i + 1, N)]
pw = lam ** (N - np.arange(N + 1))
om = np.exp(2j * np.pi / 3)
chars = om ** ((A @ A.T) % 3)          # chars[v_index, a_index] = omega^{v.a}


def gammas():
    if mode == 'all':
        for e in itertools.product(range(3), repeat=len(pairs)):
            G = np.zeros((N, N), dtype=np.int64)
            for (i, j), v in zip(pairs, e):
                G[i, j] = G[j, i] = v
            yield G
    else:
        for _ in range(count):
            G = np.zeros((N, N), dtype=np.int64)
            for (i, j) in pairs:
                G[i, j] = G[j, i] = rng.integers(0, 3)
            yield G


worst = (0.0, None); minS0 = np.inf; ties = 0
for G in gammas():
    GA = (A @ G.T) % 3                                     # Gamma a for every a
    for d in A:                                            # every shift d
        w = ((A != 0) | ((GA + d) % 3 != 0)).sum(1)
        f = pw[w]
        S = chars @ f                                      # S(d, v) for all v
        S0 = S[0].real
        minS0 = min(minS0, S0)
        ratio = np.abs(S[1:]).max() / S0
        if ratio > worst[0] + 1e-12:
            worst = (ratio, (G.copy(), d.copy()))
        ties += int((np.abs(np.abs(S[1:]) - S0) < 1e-9 * S0).sum())
print(f"N={N} t={t} ({mode}): max |S(d,v)|/S(d,0) over v != 0 = {worst[0]:.10f}; min S(d,0) = {minS0:.4f}; equality cases {ties}")
