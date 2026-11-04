"""Checks for Lemma 13 (bath reduction) and Proposition 14 (balanced triple) of the README.

p = p_G (x) |0>^{(x)m}: p_G a single excitation on the k group copies, m bath copies in the vacuum.  A = M^Gamma and
B = R_{p_G}^Gamma on the group.  Lemma 13 writes Y_c(p) as a direct sum of blocks:
- N_2 = 1 (one bath copy excited in row and column, letters (i, i')): 2^{-N_0} [A00 (x) (I - F/2) - B00 (x) F/2] on G0_j;
- N_2 = 0 (N bath copies carry one column excitation; the row excitation sits in G or on one of them):
  2^{-N_0} F_N with F_N = [[Y11, -1/2 Y10 (x) 1^T], [-1/2 Y01 (x) 1, A00 (x) I_N + 1/4 Y00 (x) (J_N - I_N)]].
(a) For random group vectors (k = 2, 3) and m = 1, 2, 3 bath copies, the multiset of eigenvalues of all blocks equals
    the spectrum of Y_c(p) computed directly by relative_anchor.py (every c).
(b) For the balanced triple |100> + |010> + |001>, Y_c(p) >= 0 directly for n = 3, ..., 6 and every c (Proposition 14;
    the exact proof is check_balanced_triple_exact.py)."""
import numpy as np, itertools
from math import comb
from relative_anchor import Relative, Mel

def group_ops(pG, k):
    """A and B = R_{pG}^Gamma on the group, as functions of ((row, col'), (row', col)) for rows with <= 1 excitation."""
    S = list(itertools.product(range(3), repeat=k)); pos = {s: i for i, s in enumerate(S)}
    P = list(pG.items())
    D = np.array([[sum(np.conj(pa) * pb * Mel(a, x, b, xp) for a, pa in P for b, pb in P) for xp in S] for x in S])
    Di = np.linalg.inv(D)
    rows = [s for s in S if sum(v != 0 for v in s) <= 1]
    MV = {(b, y): np.array([sum(Mel(b, y, a, z) * pa for a, pa in P) for z in S]) for b in rows for y in S}
    def Aent(a, yp, b, y): return Mel(b, yp, a, y)
    def Bent(a, yp, b, y): return MV[(b, yp)] @ Di @ MV[(a, y)].conj()
    return S, Aent, Bent

def block_spectrum(pG, k, m, c):
    S, Aent, Bent = group_ops(pG, k)
    q = lambda s: sum(v != 0 for v in s); vac = (0,) * k
    H = {j: [s for s in S if q(s) == j] for j in range(k + 1)}
    def mat(f, X, Y): return np.array([[f(a, yp, b, y) for (b, y) in Y] for (a, yp) in X])
    G1 = lambda j: [(a, y) for a in H[1] for y in H[j]] if 0 <= j <= k else []
    G0 = lambda j: [(vac, y) for y in H[j]] if 0 <= j <= k else []
    ev = []
    for N2 in (0, 1):
        for N1 in range(0, m - N2 + 1):
            N0 = m - N2 - N1
            mult = comb(m, N2) * comb(m - N2, N1) * 2 ** N1          # positions and letters of the charge-1 copies
            if N2 == 1:
                j = c - N1 - 1; X = G0(j)
                if not X: continue
                A00 = mat(Aent, X, X); B00 = mat(Bent, X, X)
                # letters (i, i') at the doubly excited copy: F = +1 on 3 symmetric states, -1 on 1 antisymmetric
                for sgn, cnt in ((1, 3), (-1, 1)):
                    blk = A00 * (1 - sgn / 2) - B00 * (sgn / 2)
                    ev += list(np.repeat(np.linalg.eigvalsh((blk + blk.conj().T) / 2) * 2.0 ** -N0, cnt * mult))
            else:
                j = c - N1; X1 = G1(j); X0 = G0(j + 1) if N1 >= 1 else []
                if not X1 and not X0: continue
                N = N1
                Y = lambda X, Z: mat(Aent, X, Z) - mat(Bent, X, Z)
                if X0:
                    low = np.kron(mat(Aent, X0, X0), np.eye(N)) + 0.25 * np.kron(Y(X0, X0), np.ones((N, N)) - np.eye(N))
                if X1 and X0:
                    off = -0.5 * np.kron(Y(X1, X0), np.ones((1, N)))
                    F = np.block([[Y(X1, X1), off], [off.conj().T, low]])
                elif X1:
                    F = Y(X1, X1)
                else:
                    F = low
                ev += list(np.repeat(np.linalg.eigvalsh((F + F.conj().T) / 2) * 2.0 ** -N0, mult))
    return np.sort(np.array(ev).real)

if __name__ == '__main__':
    rng = np.random.default_rng(13); ok = True
    for k in (2, 3):
        for m in (1, 2, 3):
            n = k + m
            for c in range(1, n + 1):
                w = rng.random(k) + 0.1
                pG = {}
                for t in range(k):
                    s = [0] * k; s[t] = 1; pG[tuple(s)] = w[t]
                R = Relative(n, c); p = R.pvec(np.concatenate([w, np.zeros(m)]))
                direct = np.sort(np.linalg.eigvalsh(R.Y(p)[0]))
                blocks = block_spectrum(pG, k, m, c)
                err = np.abs(direct - blocks).max() if len(direct) == len(blocks) else np.inf
                ok &= err < 1e-10
                if c in (1, n // 2 + 1, n):
                    print(f"(a) k={k} m={m} c={c}: {len(direct)} eigenvalues, max deviation {err:.1e}", flush=True)
    worst = 9
    for n in (3, 4, 5, 6):
        for c in range(1, n + 1):
            R = Relative(n, c); worst = min(worst, np.linalg.eigvalsh(R.Y(R.pvec(np.array([1.0, 1, 1] + [0] * (n - 3))))[0])[0])
    print(f"(b) balanced triple, n = 3..6, every c: min eigenvalue of Y_c {worst:+.2e}"); ok &= worst > -1e-12
    print("PASS" if ok else "FAIL")
