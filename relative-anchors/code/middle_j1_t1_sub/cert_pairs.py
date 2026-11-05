"""Exact pair-decomposition certificates for the gauge-fixed |T| = 1 block of (R_inf)_1 (NOTES.md, Section 4).

For a rational weight vector w and the copy t, take the block from the Lemma B builder j1lib.py (exact Fractions),
delete the state P_t (gauge fixing; the kernel vector p (x) |2@t> has P_t-coordinate w_t != 0) and group the states by
copies: the group of t is {R_t}, the group of c != t is {R_c, P_c, Q_c}.  A certificate is a list of rational symmetric
matrices Z_ab, one for every pair of groups (a, b), supported on the states of a and b, each positive definite, and
one positive definite matrix S_a for every group, such that

    G = sum_{a<b} Z_ab + sum_a S_a        exactly.

It proves that G is positive definite (so the block is >= 0 with kernel exactly p (x) |2@t>) at that weight vector,
and it exhibits G as a sum of terms each involving the states of at most two copies.  The certificate is found by a
floating-point semidefinite program (cvxpy/SCS: maximise eps with G - eps diag(G) in the cone), then made exact:
the off-diagonal block of Z_ab is set to the exact block of G, its diagonal blocks are rounded to rationals with a
margin, and the remainders S_a are computed exactly.  Only the exact checks (LDL with positive pivots, exact sum) are
claimed; the SDP is a search tool.

Negative control: with the B-part scaled by 1001/1000 (j1lib bscale) at points next to a degenerate corner, the
scaled gauge-fixed block is not positive definite (exact LDL), so no certificate can exist; the SDP returns eps < 0
there and the rounding step fails the exact check.

Usage: python cert_pairs.py [seed]"""
import sys
import random
import itertools
from fractions import Fraction as Fr
import numpy as np
import cvxpy as cp
from j1lib import J1, block_states


def ldl_pd(A):
    A = [row[:] for row in A]
    n = len(A)
    for c in range(n):
        if A[c][c] <= 0:
            return False
        for i in range(c + 1, n):
            if A[i][c] != 0:
                f = A[i][c] / A[c][c]
                for j in range(c + 1, n):
                    A[i][j] -= f * A[c][j]
    return True


def gauge_fixed(w, t, bscale=1):
    k = len(w)
    st = block_states(k, (t,))
    M = J1(w, bscale=bscale).matrix(st)
    def find(pred):
        return [i for i, s in enumerate(st) if pred(s)][0]
    groups = []
    for c in range(k):
        if c == t:
            groups.append([find(lambda s: s[0] == (2, t) and s[1][t] == 1)])
        else:
            groups.append([find(lambda s: s[0] == (2, t) and s[1][c] == 1),
                           find(lambda s: s[0] == (1, c)),
                           find(lambda s: s[0] is None and s[1][c] == 1)])
    pt = find(lambda s: s[0] == (1, t))
    keep = [i for i in range(len(st)) if i != pt]
    pos = {i: n for n, i in enumerate(keep)}
    G = [[M[i][j] for j in keep] for i in keep]
    groups = [[pos[i] for i in g] for g in groups]
    return G, groups


def certify(G, groups, den=10 ** 12):
    n = len(G)
    Gf = np.array([[float(v) for v in row] for row in G])
    dg = np.sqrt(np.diag(Gf))
    Gn = Gf / dg[:, None] / dg[None, :]
    pairs = list(itertools.combinations(range(len(groups)), 2))
    Zs = {}
    tot = 0
    for a, b in pairs:
        idx = groups[a] + groups[b]
        Z = cp.Variable((len(idx), len(idx)), PSD=True)
        E = np.zeros((len(idx), n))
        for i, j in enumerate(idx):
            E[i, j] = 1
        tot = tot + E.T @ Z @ E
        Zs[(a, b)] = (Z, idx)
    eps = cp.Variable()
    prob = cp.Problem(cp.Maximize(eps), [tot + eps * np.eye(n) == Gn, eps <= 1])
    prob.solve(solver=cp.SCS, eps=1e-10, max_iters=200000)
    e = float(eps.value) if eps.value is not None else -1.0
    if e <= 0:
        return False, e
    # margin: shift every Z_ab by delta on its diagonal (normalised frame), take it from the slack eps
    cnt = {a: sum(1 for p in pairs if a in p) for a in range(len(groups))}
    delta = e / (4 * max(cnt.values()))
    rat = lambda v: Fr(int(round(v * den)), den)
    Zr = {}
    diag_used = [[Fr(0)] * n for _ in range(n)]
    for (a, b), (Z, idx) in Zs.items():
        Zn = Z.value + delta * np.eye(len(idx))
        Zb = Zn * dg[idx][:, None] * dg[idx][None, :]          # back to the frame of G
        m = len(idx)
        R = [[None] * m for _ in range(m)]
        la = len(groups[a])
        for i in range(m):
            for j in range(m):
                same = (i < la) == (j < la)
                if same:
                    if j >= i:
                        R[i][j] = rat(Zb[i, j])
                    else:
                        R[i][j] = R[j][i]
                else:
                    R[i][j] = G[idx[i]][idx[j]]                 # exact off-diagonal block
        Zr[(a, b)] = (R, idx)
        for i in range(m):
            for j in range(m):
                if (i < la) == (j < la):
                    diag_used[idx[i]][idx[j]] += R[i][j]
    ok = all(ldl_pd(R) for R, idx in Zr.values())
    S = {}
    for a, g in enumerate(groups):
        Sa = [[G[i][j] - diag_used[i][j] for j in g] for i in g]
        S[a] = Sa
        ok &= ldl_pd(Sa)
    # exact reconstruction
    tot = [[Fr(0)] * n for _ in range(n)]
    for R, idx in Zr.values():
        for i, ii in enumerate(idx):
            for j, jj in enumerate(idx):
                tot[ii][jj] += R[i][j]
    for a, g in enumerate(groups):
        for i, ii in enumerate(g):
            for j, jj in enumerate(g):
                tot[ii][jj] += S[a][i][j]
    ok &= tot == G
    return ok, e


if __name__ == '__main__':
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    random.seed(seed)
    allok = True
    pts = []
    for k in list(range(3, 9)) + [10]:
        for trial in range(4 if k <= 5 else 2):
            mode = trial % 4
            if mode == 0:
                w = [Fr(random.randint(1, 20), random.randint(1, 20)) for _ in range(k)]
            elif mode == 1:
                w = [Fr(1, 3 ** random.randint(0, 4)) for _ in range(k)]
            elif mode == 2:
                w = [Fr(random.randint(5, 20), 10) for _ in range(k)]
                w[random.randrange(k)] = Fr(1, 100)
            else:
                w = [Fr(random.randint(1, 4), 100) for _ in range(k)]
                w[random.randrange(k)] = Fr(1)
            pts.append((w, random.randrange(k)))
    # the two hard limits (Section 5 of middle_j1_t1/NOTES.md): small weight on t with two dominant U-copies,
    # and t carrying almost all the weight
    pts.append(([Fr(1, 100), Fr(1), Fr(1)], 0))
    pts.append(([Fr(1, 60), Fr(1), Fr(1), Fr(1, 20)], 0))
    pts.append(([Fr(30), Fr(1), Fr(1, 7)], 0))
    pts.append(([Fr(40), Fr(1), Fr(1, 3), Fr(1, 3)], 0))
    # the shape of the floating-point near-failure of num_pairs.py (Section 4), with x_t about 1e-4
    pts.append(([Fr(1, 100), Fr(36, 100), Fr(1, 300), Fr(17, 1000), Fr(18, 100), Fr(17, 100), Fr(22, 1000), Fr(90, 100),
                 Fr(35, 1000)], 0))
    for w, t in pts:
        G, groups = gauge_fixed(w, t)
        ok, e = certify(G, groups)
        allok &= ok
        print(f'k = {len(w)}, t = {t}, w = {[str(a) for a in w]}: SDP margin {e:.3e}; exact certificate '
              f'{"verified" if ok else "FAILED"}', flush=True)
    print(f'pair certificates at {len(pts)} rational weight vectors: {"PASS" if allok else "FAIL"}', flush=True)
    # negative controls
    negok = True
    for w, t in [([Fr(1, 300), Fr(1), Fr(1)], 0), ([Fr(1, 200), Fr(1), Fr(1), Fr(1, 3)], 0),
                 ([Fr(1), Fr(2), Fr(1, 400), Fr(1), Fr(1, 2)], 2)]:
        G1, gr = gauge_fixed(w, t)
        G2, gr2 = gauge_fixed(w, t, bscale=Fr(1001, 1000))
        pd1, pd2 = ldl_pd(G1), ldl_pd(G2)
        ok2, e2 = certify(G2, gr2)
        print(f'  control w = {[str(a) for a in w]}, t = {t}: true block PD {pd1}; scaled block PD {pd2}, SDP margin '
              f'{e2:.3e}, certificate {"verified" if ok2 else "fails"}', flush=True)
        negok &= pd1 and (not pd2) and (not ok2)
    print(f'negative controls (B scaled by 1001/1000): {"PASS" if negok else "FAIL"}', flush=True)
    print('ALL PASS' if (allok and negok) else 'SOME FAIL')
