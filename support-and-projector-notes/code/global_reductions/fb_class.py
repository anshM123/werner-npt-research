"""Merging budget FB(t) on the stabilizer pairs C = a C1 + b C2 (Section 7 class), sites i, j merged.
Q = S(d,0)/3^N, X = S(d,v)/3^N (all N sites); Q', X' the same sums with the merged weight (1 - 9t if the Pauli is the
identity on both merged sites, else 1) times lam^{#identity sites among the others}.
FB on the class <=> for all d, v != 0:  [Q - g Q' - (1-t)^{N-1}/2] - [|X - g X'| - t(1-t)^{N-2}/2] >= 0  (|a| = |b|)
and Q - g Q' - (1-t)^{N-1}/2 >= 0 (rank one), with g = (1-2t)/2.  Reports the minima of both, divided by Q.
usage: python fb_class.py N t family i j"""
import sys, itertools
import numpy as np
N, t, fam, mi, mj = int(sys.argv[1]), float(sys.argv[2]), sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
lam = 1 - 3 * t; gam = (1 - 2 * t) / 2
G = np.zeros((N, N), dtype=np.int64)
def edge(a, b): G[a, b] = G[b, a] = 1
if fam == 'star':
    for i in range(1, N): edge(0, i)
elif fam == 'path':
    for i in range(N - 1): edge(i, i + 1)
elif fam == 'cycle':
    for i in range(N): edge(i, (i + 1) % N)
elif fam.startswith('bip'):
    h = int(fam[3:])
    for i in range(h):
        for j in range(h, N): edge(i, j)
elif fam.startswith('rand'):
    rng = np.random.default_rng(int(fam[4:]))
    M = np.triu(rng.integers(0, 3, size=(N, N)), 1); G = (M + M.T) % 3
A = np.array(list(itertools.product(range(3), repeat=N)), dtype=np.int8)
GA = (A.astype(np.int64) @ G.T) % 3; nzA = A != 0
others = [s for s in range(N) if s not in (mi, mj)]
worst2, worst1 = np.inf, np.inf; arg = None
for d in A:
    idn = ~(nzA | ((GA + d) % 3 != 0))                         # identity at site s
    k_all = idn.sum(1)
    f = lam ** k_all
    both = idn[:, mi] & idn[:, mj]
    k_oth = idn[:, others].sum(1)
    fm = np.where(both, 1 - 9 * t, 1.0) * lam ** k_oth
    S = np.fft.fftn(f.reshape((3,) * N)).ravel() / 3 ** N
    Sm = np.fft.fftn(fm.reshape((3,) * N)).ravel() / 3 ** N
    Q, Qm = S[0].real, Sm[0].real
    diag = Q - gam * Qm - (1 - t) ** (N - 1) / 2
    cross = np.abs(S[1:] - gam * Sm[1:]).max() - t * (1 - t) ** (N - 2) / 2
    worst1 = min(worst1, diag / Q)
    if (diag - cross) / Q < worst2:
        worst2 = (diag - cross) / Q; arg = d.copy()
print(f"N={N} t={t} {fam} merge({mi},{mj}): min rank-one term {worst1:+.5f}, min pair margin {worst2:+.5f} at d={arg.tolist()}")
