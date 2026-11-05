"""Exact checks of the level reductions of lev.py against the Lemma B builder j1lib.py.

For l = 1 (equal U-weights), l = 2 (two U-weight levels; also with u_1 = 1 or u_1 = u_2 = 1 fixed), l = 3 (two
single U-copies and u_3 equal ones) and l = 4 with four single copies (k = 5, general weights), at rational weights and
level sizes u_i, with t = 0:
 (a) the P-first reduction: the coupling B of the P-states to (R, Q) satisfies Pi diag(1/N) B = B and w^T B = 0,
     so with K = [[I, -diag(1/N) B], [0, I]] (order P, then R and Q) K^T M K = diag(Pi, Phi), Phi = C - B^T diag(1/N) B;
 (b) the invariant and standard parts of D Phi D (divided by L = |s|) equal Finv and Fstd of lev.py at that point;
 (c) the same for the gauge-fixed block (state P_t deleted) and Ginv, Gstd, with the B-part scaled by 1 and by
     1 + 1/1000 (j1lib bscale);
 (d) Phi is positive definite at the point (exact LDL).
Usage: python check_lev.py"""
import sys, random, time
import flint
from fractions import Fraction as Fr
from j1lib import J1, block_states, mu
from lev import build

random.seed(11)
t0 = time.time()


def ev(f, vals):
    num = sum((Fr(int(flint.fmpq(c).p), int(flint.fmpq(c).q)) * prod(vals, m) for m, c in f.n.to_dict().items()), Fr(0))
    den = sum((Fr(int(flint.fmpq(c).p), int(flint.fmpq(c).q)) * prod(vals, m) for m, c in f.d.to_dict().items()), Fr(0))
    return num / den


def prod(vals, m):
    out = Fr(1)
    for v, e in zip(vals, m):
        if e:
            out *= v ** int(e)
    return out


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


def project(Mt, groups, levels):
    """groups: list of (name, dict level -> list of indices, t-index or None).  Invariant form on
    (group, 't' or level) coordinates; standard 2 x 2/3 x 3 forms per level."""
    coords = []
    for name, bylev, tidx in groups:
        if tidx is not None:
            coords.append([tidx])
        for i in range(levels):
            coords.append(bylev[i])
    inv = [[sum(Mt[a][b] for a in A for b in B) for B in coords] for A in coords]
    std = []
    for i in range(levels):
        S = []
        for name, bylev, tidx in groups:
            row = []
            for name2, bylev2, tidx2 in groups:
                h, hp = 0, 1
                if len(bylev[i]) < 2:
                    row = None
                    break
                row.append(Mt[bylev[i][h]][bylev2[i][h]] - Mt[bylev[i][h]][bylev2[i][hp]])
            if row is None:
                S = None
                break
            S.append(row)
        std.append(S)
    return inv, std


ok = True
cnt = 0
cases = ([(1, [2], us, None) for us in ([2], [3], [4], [6])]
         + [(2, [1, 1], us, None) for us in ([1, 2], [2, 2], [2, 3], [3, 1], [1, 4])]
         + [(2, [1, 2], us, [True, False]) for us in ([1, 2], [1, 3])]
         + [(2, [1, 1], [1, 1], [True, True])]
         + [(3, [1, 1, 2], us, [True, True, False]) for us in ([1, 1, 2], [1, 1, 3], [1, 1, 5])]
         + [(4, [1, 1, 1, 1], [1, 1, 1, 1], [True] * 4)])
SYM = {}
for l, offs, us, fixed in cases:
    key = (l, tuple(offs), None if fixed is None else tuple(fixed))
    if key not in SYM:
        SYM[key] = {bsc: build(l, offs, bsc, want_phi=(bsc == 1), fixed=fixed) for bsc in (1, Fr(1001, 1000))}
    for trial in range(2):
        wt = Fr(random.randint(1, 9), random.randint(1, 9))
        wl = [Fr(random.randint(1, 9), random.randint(1, 9)) for _ in range(l)]
        w = [wt] + [wl[i] for i in range(l) for _ in range(us[i])]
        k = len(w)
        lev_of = {}
        c = 1
        for i in range(l):
            for _ in range(us[i]):
                lev_of[c] = i
                c += 1
        p = wt * wt
        qs = [v * v for v in wl]
        L = p + sum(us[i] * qs[i] for i in range(l))
        vals = [p] + qs + [Fr(us[i] - offs[i]) for i in range(l) if not (fixed and fixed[i])]
        st = block_states(k, (0,))
        R = [i for i, (a, y) in enumerate(st) if a is not None and a[0] == 2]
        P = [i for i, (a, y) in enumerate(st) if a is not None and a[0] == 1]
        Q = [i for i, (a, y) in enumerate(st) if a is None]
        for bsc in (1, Fr(1001, 1000)):
            J = J1(w, bscale=bsc)
            M = J.matrix(st)
            m = [mu(z, J.w) for z in st]
            sym = SYM[key][bsc]
            if bsc == 1:
                # (a) P-first reduction
                N = [Fr(1) if r == 0 else Fr(2) for r in range(k)]
                E = R + Q
                Bm = [[M[P[r]][e] for e in E] for r in range(k)]
                Pi = [[M[P[r]][P[r2]] for r2 in range(k)] for r in range(k)]
                okA = all(sum(Pi[r][r2] * Bm[r2][e] / N[r2] for r2 in range(k)) == Bm[r][e] for r in range(k) for e in range(len(E)))
                okA &= all(sum(w[r] * Bm[r][e] for r in range(k)) == 0 for e in range(len(E)))
                Phi = [[M[E[a]][E[b2]] - sum(Bm[r][a] * Bm[r][b2] / N[r] for r in range(k)) for b2 in range(len(E))] for a in range(len(E))]
                mE = [m[e] for e in E]
                Pt = [[mE[a] * Phi[a][b2] * mE[b2] / L for b2 in range(len(E))] for a in range(len(E))]
                # groups in Phi coordinates: R (t + levels), Q (levels)
                Rl = {i: [a for a in range(k) if a != 0 and lev_of[a] == i] for i in range(l)}
                Ql = {i: [k + j for j, e in enumerate(Q) if lev_of[st[e][1].index(1)] == i] for i in range(l)}
                inv, std = project(Pt, [('R', Rl, 0), ('Q', Ql, None)], l)
                okB = all(inv[a][b2] == ev(sym['Finv'][a][b2], vals) for a in range(len(inv)) for b2 in range(len(inv)))
                for i in range(l):
                    if std[i] is not None:
                        okB &= all(std[i][a][b2] == ev(sym['Fstd'][i][a][b2], vals) for a in range(2) for b2 in range(2))
                okD = ldl_pd(Phi)
            # (c) gauge-fixed block, P_t deleted
            G = R + [P[r] for r in range(1, k)] + Q
            Mt = [[m[a] * M[a][b2] * m[b2] / L for b2 in G] for a in G]
            Rl = {i: [a for a in range(k) if a != 0 and lev_of[a] == i] for i in range(l)}
            Pl = {i: [k + r - 1 for r in range(1, k) if lev_of[r] == i] for i in range(l)}
            Ql = {i: [2 * k - 1 + j for j, e in enumerate(Q) if lev_of[st[e][1].index(1)] == i] for i in range(l)}
            inv, std = project(Mt, [('R', Rl, 0), ('P', Pl, None), ('Q', Ql, None)], l)
            okC = all(inv[a][b2] == ev(sym['Ginv'][a][b2], vals) for a in range(len(inv)) for b2 in range(len(inv)))
            for i in range(l):
                if std[i] is not None:
                    okC &= all(std[i][a][b2] == ev(sym['Gstd'][i][a][b2], vals) for a in range(3) for b2 in range(3))
            good = okC and (okA and okB and okD if bsc == 1 else True)
            if trial == 0:
                print(f'   bscale {bsc}: invariant gauge form size {len(inv)}, std forms {[None if S is None else len(S) for S in std]}, '
                      f'G_inv[0][0] = {inv[0][0]} = {ev(sym["Ginv"][0][0], vals)}', flush=True)
            if not good:
                print(f'FAIL l={l} us={us} w={[str(v) for v in w]} bscale={bsc}: '
                      f'{(okA, okB, okD) if bsc == 1 else ""} {okC}', flush=True)
            ok &= good
            cnt += 1
    print(f'l = {l}, u = {us}, fixed {fixed}: done, all pass so far: {ok} [{time.time() - t0:.1f}s]', flush=True)
print(f'{cnt} exact comparisons')
print('PASS' if ok else 'FAIL')
