"""Check of Lemma 13' (README): for p = p_G (x) |0>^{(x)m}, the tangent compressions of Lemma 11' (rows H_0 + H_1, all
columns) split into blocks whose congruence classes are Y|_{G1_j + G0_{j+1}} (+) (1/4)(3A + B)^{00} and the blocks of
Lemma 13; so they are >= 0 for all m iff (R_inf) and (R2) hold.  Here we compare the multiset of eigenvalues of the
actual blocks with the direct spectrum of (M - R_p)^Gamma restricted to (H_0 + H_1) (x) (all columns)."""
import numpy as np, itertools
from math import comb
from relative_anchor import Mel
from check_bath_reduction import group_ops

def direct_spectrum(w, m):
    k = len(w); n = k + m
    S = list(itertools.product(range(3), repeat=n)); pos = {s: i for i, s in enumerate(S)}
    H1 = [s for s in S if sum(x != 0 for x in s) == 1]; vac = (0,) * n
    p = {}
    for t in range(k):
        s = [0] * n; s[t] = 1; p[tuple(s)] = w[t]
    P = list(p.items())
    D = np.array([[sum(pa * pb * Mel(a, x, b, xp) for a, pa in P for b, pb in P) for xp in S] for x in S])
    Di = np.linalg.inv(D)
    rows = [vac] + H1
    idx = [(a, y) for a in rows for y in S]; ui = {q: i for i, q in enumerate(idx)}
    G = np.array([[sum(Mel(a, y, b, x) * pb for b, pb in P) for x in S] for (a, y) in idx])
    Mc = np.array([[Mel(a, y, b, yp) for (b, yp) in idx] for (a, y) in idx])
    W = Mc - G @ Di @ G.T
    P1 = np.array([[ui[(b, yp)] for (b, y) in idx] for (a, yp) in idx]); P2 = np.array([[ui[(a, y)] for (b, y) in idx] for (a, yp) in idx])
    Y = W[P1, P2]; return np.sort(np.linalg.eigvalsh((Y + Y.T) / 2))

def block_spectrum(w, m):
    k = len(w); pG = {}
    for t in range(k):
        s = [0] * k; s[t] = 1; pG[tuple(s)] = w[t]
    S, Aent, Bent = group_ops(pG, k)
    q = lambda s: sum(v != 0 for v in s); vac = (0,) * k
    H = {j: [s for s in S if q(s) == j] for j in range(k + 1)}
    G1 = lambda j: [(a, y) for a in H[1] for y in H[j]] if 0 <= j <= k else []
    G0 = lambda j: [(vac, y) for y in H[j]] if 0 <= j <= k else []
    mat = lambda f, X, Z: np.array([[f(a, yp, b, y) for (b, y) in Z] for (a, yp) in X]).reshape(len(X), len(Z))
    Y = lambda X, Z: mat(Aent, X, Z) - mat(Bent, X, Z)
    ev = []
    for N2 in (0, 1):
        for N in range(0, m - N2 + 1):
            N0 = m - N2 - N; mult = comb(m, N2) * comb(m - N2, N) * 2 ** N
            if N2 == 1:
                for j in range(0, k + 1):
                    X = G0(j); A00 = mat(Aent, X, X); B00 = mat(Bent, X, X)
                    ev += list(np.repeat(np.linalg.eigvalsh(0.5 * (A00 - B00)), 3 * mult) * 2.0 ** -N0)
                    ev += list(np.repeat(np.linalg.eigvalsh(1.5 * A00 + 0.5 * B00), mult) * 2.0 ** -N0)
                continue
            for j in range(-1, k + 1):
                X1, X0 = G1(j), G0(j + 1)
                if not X1 and not X0: continue
                if N == 0:
                    F = Y(X1 + X0, X1 + X0)
                else:
                    a, b = len(X1), len(X0)
                    dim = a + b * N + b; F = np.zeros((dim, dim))
                    i1, i2, i3 = slice(0, a), slice(a, a + b * N), slice(a + b * N, dim)
                    one = np.ones((1, N))
                    if a:
                        F[i1, i1] = Y(X1, X1); F[i1, i3] = Y(X1, X0); F[i3, i1] = F[i1, i3].T
                        F[i1, i2] = -0.5 * np.kron(Y(X1, X0), one); F[i2, i1] = F[i1, i2].T
                    F[i2, i2] = np.kron(mat(Aent, X0, X0), np.eye(N)) + 0.25 * np.kron(Y(X0, X0), np.ones((N, N)) - np.eye(N))
                    F[i3, i3] = Y(X0, X0)
                    F[i2, i3] = -0.5 * np.kron(Y(X0, X0), one.T); F[i3, i2] = F[i2, i3].T
                ev += list(np.repeat(np.linalg.eigvalsh((F + F.T) / 2) * 2.0 ** -N0, mult))
    return np.sort(np.array(ev))

if __name__ == '__main__':
    rng = np.random.default_rng(21); ok = True
    for k, m in [(2, 1), (2, 2), (3, 1), (3, 2), (2, 3)]:
        w = rng.random(k) + 0.1
        d = direct_spectrum(w, m); b = block_spectrum(w, m)
        err = np.abs(d - b).max() if len(d) == len(b) else np.inf
        print(f"k={k} m={m}: {len(d)} eigenvalues, max deviation {err:.1e}, smallest {d[0]:+.2e}", flush=True)
        ok &= err < 1e-10
    print("PASS" if ok else "FAIL")
