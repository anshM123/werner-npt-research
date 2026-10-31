"""Fast stabiliser-code block helpers (vectorised): adapted symplectic bases, exact weights, Lagrangian subspaces, and
   a floating-point LP test for k up to 3. Exact certificates are in check_code_blocks_margin.py.
   Paulis = vectors in F_3^(2b). Choose a symplectic basis adapted to the code: stabiliser generators g_l, destabilisers
   d_l (<d_l, g_m> = delta), logical pairs (x_j, z_j). Coordinates of a Pauli: syndrome alpha (from destabilisers),
   logical beta in F_3^(2k), stabiliser part gamma. Weight p(alpha, beta) = |S|^-1 sum_gamma (-1/2)^(b - wt).
   Separable vertices: uniform 3^-k on cosets beta0 + Lag, Lag a Lagrangian subspace of F_3^(2k) (standard form)."""
import itertools, sys, numpy as np
from fractions import Fraction
from scipy.optimize import linprog
from scipy.sparse import csc_matrix


def J(n):
    return np.block([[np.zeros((n, n), int), np.eye(n, dtype=int)], [-np.eye(n, dtype=int), np.zeros((n, n), int)]]) % 3


def sform(u, v, n):
    return int(u @ J(n) @ v) % 3


def inv3(x):
    return {1: 1, 2: 2}[x % 3]


def adapted_basis(gens, b):
    """symplectic basis: list of (destab, stab) pairs and logical (x, z) pairs, over F_3."""
    Jb = J(b)
    def f(u, v): return int(u @ Jb @ v) % 3
    stabs = [np.array(g) % 3 for g in gens]
    # complete to a basis of F_3^(2b): we find vectors with prescribed pairings by solving linear systems
    def solve(rows, rhs):
        # find x with rows @ Jb @ x = rhs (mod 3), brute force via Gaussian elimination
        A = np.array([r @ Jb % 3 for r in rows]) % 3 if rows else np.zeros((0, 2 * b), int)
        M = np.concatenate([A, np.array(rhs).reshape(-1, 1)], 1) % 3 if rows else None
        x = np.zeros(2 * b, int)
        if not rows: return x
        M = M.copy(); r = 0; piv = []
        for c in range(2 * b):
            p = next((i for i in range(r, len(M)) if M[i, c] % 3), None)
            if p is None: continue
            M[[r, p]] = M[[p, r]]; M[r] = (M[r] * inv3(M[r, c])) % 3
            for i in range(len(M)):
                if i != r and M[i, c]: M[i] = (M[i] - M[i, c] * M[r]) % 3
            piv.append(c); r += 1
        assert all(M[i, -1] == 0 for i in range(r, len(M))), "inconsistent"
        for i, c in enumerate(piv): x[c] = M[i, -1]
        return x
    destabs = []
    for l in range(len(stabs)):
        rows = stabs + destabs
        rhs = [1 if m == l else 0 for m in range(len(stabs))] + [0] * len(destabs)
        destabs.append(solve(rows, rhs))
    # logical pairs: vectors orthogonal to all stabs and destabs, made symplectic
    logx, logz = [], []
    k = b - len(stabs)
    space = [np.array(v) for v in itertools.product(range(3), repeat=2 * b)]
    used = stabs + destabs
    def ok(v, lst): return all(f(v, w) == 0 for w in lst)
    for _ in range(k):
        cons = stabs + destabs + logx + logz
        x = next(v for v in space if any(v) and ok(v, cons) and
                 np.linalg.matrix_rank(np.array(cons + [v], float)) > len(cons))
        # z: pairing 1 with x, 0 with the rest
        z = solve(cons + [x], [0] * len(cons) + [1])
        logx.append(x); logz.append(z)
    return stabs, destabs, logx, logz


def lagrangians(k):
    """all Lagrangian subspaces of F_3^(2k) with the standard form, as sorted tuples of codes."""
    vecs = np.array(list(itertools.product(range(3), repeat=2 * k)))
    code = vecs @ (3 ** np.arange(2 * k)[::-1])
    Jk = J(k)
    G = (vecs @ Jk @ vecs.T) % 3
    level = {(0,)}
    for d in range(k):
        nxt = set()
        for sub in level:
            sub_idx = list(sub)
            okv = np.all(G[sub_idx][:, :] == 0, axis=0)
            for v in np.nonzero(okv)[0]:
                if v in sub: continue
                span = {int(((vecs[s] + a * vecs[v]) % 3) @ (3 ** np.arange(2 * k)[::-1])) for s in sub for a in range(3)}
                nxt.add(tuple(sorted(span)))
        level = nxt
    return [list(L) for L in level], vecs


def run(name, gens, b, lp=True, exact=False):
    stabs, destabs, logx, logz = adapted_basis(gens, b)
    k = len(logx)
    B = np.array(destabs + logx + logz + stabs) % 3          # rows: alpha | beta(x) | beta(z) | gamma coordinates
    # coordinates: P = sum alpha_l d_l + sum beta_j (x_j or z_j) + sum gamma_l g_l; invert over F_3 by brute force
    m = len(stabs)
    allP = np.array(list(itertools.product(range(3), repeat=2 * b)))
    wt = ((allP[:, :b] != 0) | (allP[:, b:] != 0)).sum(1)
    # coefficients via pairing: alpha_l = <P, g_l> (since <d_l, g_m> = delta and others 0) up to sign conventions;
    Jb = J(b)
    alpha = (allP @ Jb @ np.array(stabs).T) % 3                       # <P, g_l>, determined by destabiliser part
    bx = (allP @ Jb @ np.array(logz).T) % 3 if k else np.zeros((len(allP), 0), int)   # <P, z_j> = coefficient of x_j
    bz = (-(allP @ Jb @ np.array(logx).T)) % 3 if k else np.zeros((len(allP), 0), int)  # <P, x_j> = -coef of z_j
    beta = np.concatenate([bx, bz], 1)
    acode = alpha @ (3 ** np.arange(m)[::-1]); bcode = beta @ (3 ** np.arange(2 * k)[::-1]) if k else np.zeros(len(allP), int)
    # weight sums
    val = np.array([Fraction(-1, 2) ** int(b - w) for w in range(b + 1)], dtype=object)
    p = {}
    for a, bb, w in zip(acode, bcode, wt):
        p[(a, bb)] = p.get((a, bb), Fraction(0)) + val[w]
    S = 3 ** m
    p = {key: v / S for key, v in p.items()}
    assert all(len([1 for key in p if key[0] == a]) == 9 ** k for a in set(acode))
    vals = sorted(set(p.values()))
    print(f"{name}: b = {b}, k = {k}, |S| = {S}; weights {[str(v) for v in vals]}; min {vals[0]}")
    if not lp: return
    lags, vecs = lagrangians(k)
    print(f"  {len(lags)} Lagrangian subspaces of F_3^{2 * k}")
    nb = 9 ** k
    # cosets: for each Lagrangian and each translate
    cols = set()
    pw = 3 ** np.arange(2 * k)[::-1]
    for L in lags:
        Lv = vecs[L]
        for t in range(nb):
            cs = tuple(sorted(int(c) for c in ((Lv + vecs[t]) % 3) @ pw))
            cols.add(cs)
    cols = sorted(cols)
    rows, cc = [], []
    for j, col in enumerate(cols):
        rows.extend(col); cc.extend([j] * len(col))
    A = csc_matrix((np.full(len(rows), 1.0 / 3 ** k), (rows, cc)), shape=(nb, len(cols)))
    allok = True
    for a in sorted(set(acode)):
        target = np.array([float(p[(a, t)]) for t in range(nb)])
        res = linprog(np.zeros(len(cols)), A_eq=A, b_eq=target, bounds=(0, None), method="highs")
        ok = res.status == 0; allok &= ok
        print(f"  sector {a}: p in [{min(p[(a, t)] for t in range(nb))}, {max(p[(a, t)] for t in range(nb))}], "
              f"{len(cols)} cosets, LP {'feasible' if ok else 'INFEASIBLE: ' + res.message}")
    print(f"  separable by stabiliser cosets in every sector: {allok}")


if __name__ == "__main__":
    def Xall(n): v = np.zeros(2 * n, int); v[:n] = 1; return v
    def Zc(n, c): v = np.zeros(2 * n, int); v[n:] = c; return v
    run("[[4,2,2]]", [Xall(4), Zc(4, [1, 2, 1, 2])], 4, lp=False)
    g = [np.array([1, 0, 0, 2, 0, 0, 1, 2, 0, 0])]
    for _ in range(3): g.append(np.concatenate([np.roll(g[-1][:5], 1), np.roll(g[-1][5:], 1)]))
    run("[[5,1,3]]", g, 5, lp=False)
    run("[[5,3,2]]", [Xall(5), Zc(5, [1, 1, 1, 1, 2])], 5, lp=False)
    if len(sys.argv) > 1:
        run("[[6,4,2]]", [Xall(6), Zc(6, [1, 1, 1, 1, 1, 1])], 6, lp=False)
