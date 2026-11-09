"""Cut-rank form for planes spanned by a qutrit graph state |g> and Z^u|g>:
   q_N(1/2, Pi)/2 = sum_J (-1/2)^{|J|} 3^{-r_J} (1 + chi_J(u)),  r_J = rank_F3 Gamma[J, J^c],
   chi_J(u) = [u restricted to J^c lies in the row space of Gamma[J, J^c]].
Validate against the a-sum (3^N terms) for small N, then anneal over (Gamma, u) for larger N.
usage: python cutrank.py mode N steps restarts seed   (mode = check or search)"""
import sys, itertools
import numpy as np

mode, N, steps, restarts, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
rng = np.random.default_rng(seed)
subsets = [tuple(i for i in range(N) if (m >> i) & 1) for m in range(2 ** N)]


def rref3(M):
    M = M.copy() % 3; r = 0; piv = []
    for c in range(M.shape[1]):
        p = next((i for i in range(r, M.shape[0]) if M[i, c] % 3), None)
        if p is None:
            continue
        M[[r, p]] = M[[p, r]]
        M[r] = (M[r] * (1 if M[r, c] == 1 else 2)) % 3
        for i in range(M.shape[0]):
            if i != r and M[i, c]:
                M[i] = (M[i] - M[i, c] * M[r]) % 3
        piv.append(c); r += 1
    return M[:r], piv


def in_rowspace(R, piv, v):
    v = v.copy() % 3
    for row, c in zip(R, piv):
        if v[c]:
            v = (v - v[c] * row) % 3
    return not v.any()


def value(Gam, u):
    tot = 0.0
    for J in subsets:
        Jc = [i for i in range(N) if i not in J]
        if len(J) == 0:
            r, chi = 0, float(not u.any())
        elif len(Jc) == 0:
            r, chi = 0, 1.0
        else:
            R, piv = rref3(Gam[np.ix_(J, Jc)])
            r = len(piv)
            chi = float(in_rowspace(R, piv, u[Jc]))
        tot += (-0.5) ** len(J) * 3.0 ** (-r) * (1 + chi)
    return tot          # = q_N(1/2, Pi) / 2 = q / ||Pi||^2


def asum(Gam, u):
    A = np.array(list(itertools.product(range(3), repeat=N)), dtype=np.int64)
    Ga = (A @ Gam.T) % 3
    w = ((A != 0) | (Ga != 0)).sum(1)
    inS = (A @ u) % 3 == 0
    vals = (-0.5) ** (N - w)
    return (4 * vals[inS].sum() + vals[~inS].sum()) / (2 * 3 ** N)


def rand_gamma():
    M = np.triu(rng.integers(0, 3, size=(N, N)), 1)
    return (M + M.T) % 3


if mode == 'check':
    worst = 0
    for trial in range(30):
        Gam = rand_gamma(); u = rng.integers(0, 3, size=N)
        worst = max(worst, abs(value(Gam, u) - asum(Gam, u)))
    print(f"N={N}: cut-rank form vs a-sum, max difference {worst:.2e}")
else:
    overall = np.inf
    for r in range(restarts):
        Gam = rand_gamma(); u = rng.integers(0, 3, size=N)
        cur = value(Gam, u); best = cur
        for it in range(steps):
            T = 0.02 * (1 - it / steps) + 1e-6
            Gam2, u2 = Gam.copy(), u.copy()
            if rng.random() < 0.75:
                i, j = rng.choice(N, size=2, replace=False); v = rng.integers(0, 3); Gam2[i, j] = v; Gam2[j, i] = v
            else:
                u2[rng.integers(0, N)] = rng.integers(0, 3)
            v2 = value(Gam2, u2)
            if v2 < cur or rng.random() < np.exp(-(v2 - cur) / T):
                Gam, u, cur = Gam2, u2, v2
                best = min(best, cur)
        overall = min(overall, best)
        print(f"N={N} restart {r}: min q/||Pi||^2 = {best:+.6e}", flush=True)
    print(f"N={N}: overall min over graph-state planes = {overall:+.6e}")
