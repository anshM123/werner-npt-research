"""Exact checks of the family-1 closed forms and reductions (t2lev.py) against the Lemma B builder j1lib.py (j = 2).

(a) Closed forms, entry by entry: with every copy of U in its own level (levels of size 1), the reduced matrices of
    the invariant-operator algebra are the full rescaled blocks.  At general rational weights (k = 4, 5, 6, all
    weights distinct) every block (R-R, R-R', R'-R', R-P, R'-P, P-P, R-Q, R'-Q, P-Q, Q-Q) equals the rescaled
    j1lib.py block exactly, with the B-part scaled by 1 and by 1001/1000.
(b) Lemma P2: the kernel vector p (x) |2@t 2@tau> (w on the P-states) is in the kernel, Pi w = 0, B^T w = 0 and
    Pi diag(1/N) B = B, at the same points.
(c) Two levels: at rational weights and level sizes (u0, u1), the invariant and standard reduced forms of Phi and of the
    gauge-fixed block G (P_t deleted; B scaled by 1 and 1001/1000) equal the projections of the j1lib.py block on the
    explicit test vectors (orbit indicators; e_h - e_h' on a level), exactly.
Usage: python check_t2.py"""
import itertools, random, time
from fractions import Fraction as Fr
from j1lib import J1
from r2lib import states_T2
import t2lev

t0 = time.time()
random.seed(11)
ok = True


def jblock(w, bscale=1):
    k = len(w)
    st = states_T2(k, 0, 1)
    J = J1(w, bscale=bscale)
    M = J.matrix(st)
    mu = []
    for (a, y) in st:
        if a is not None and a[0] == 1:
            mu.append(w[a[1]])
        else:
            mu.append(w[[c for c in range(k) if y[c] == 1][0]])
    Mt = [[M[i][j] * mu[i] * mu[j] for j in range(len(st))] for i in range(len(st))]
    return st, M, Mt


def index_sets(st, k):
    R = [i for i, (a, y) in enumerate(st) if a == (2, 0)]
    Rp = [i for i, (a, y) in enumerate(st) if a == (2, 1)]
    P = [i for i, (a, y) in enumerate(st) if a is not None and a[0] == 1]
    Q = [i for i, (a, y) in enumerate(st) if a is None]
    return R, Rp, P, Q


def copy_of(state):
    a, y = state
    if a is not None and a[0] == 1:
        return a[1]
    return [c for c in range(len(y)) if y[c] == 1][0]


# ---------------- (a), (b): all copies distinct
for k in (4, 5, 6):
    for trial in range(2):
        w = [Fr(random.randint(1, 9), random.randint(1, 9)) for _ in range(k)]
        s = [x * x for x in w]
        for bsc in (Fr(1), Fr(1001, 1000)):
            st, M, Mt = jblock(w, bsc)
            R, Rp, P, Q = index_sets(st, k)
            alg = t2lev.point(s[0], s[1], s[2:], [1] * (k - 2), bsc)
            B = alg['blocks']
            orb = {0: 't', 1: 'T'}
            for c in range(2, k):
                orb[c] = c - 2
            groups = dict(R=R, Rp=Rp, P=P, Q=Q)
            good = True
            for name, (gi, gj) in dict(RR=('R', 'R'), RpRp=('Rp', 'Rp'), RRp=('R', 'Rp'), PP=('P', 'P'), RP=('R', 'P'),
                                       RpP=('Rp', 'P'), QQ=('Q', 'Q'), RQ=('R', 'Q'), RpQ=('Rp', 'Q'),
                                       PQ=('P', 'Q')).items():
                op = B[name]
                for i in groups[gi]:
                    for j in groups[gj]:
                        oy, ox = orb[copy_of(st[i])], orb[copy_of(st[j])]
                        val = op.M[op.cod.index(oy)][op.dom.index(ox)]
                        if val != Mt[i][j]:
                            good = False
            ok &= good
            # (b) Lemma P2 identities (true block only)
            lp = True
            if bsc == 1:
                n = len(st)
                v = [w[st[i][0][1]] if i in P else 0 for i in range(n)]
                lp &= all(sum(M[i][j] * v[j] for j in range(n)) == 0 for i in range(n))
                Nr = {r: (1 if r in (0, 1) else 2) for r in range(k)}
                others = R + Rp + Q
                for j in others:
                    # (Pi diag(1/N) B)[r, j] = B[r, j]
                    for i in P:
                        lhs = sum(M[i][l] * M[l][j] / Nr[st[l][0][1]] for l in P)
                        lp &= lhs == M[i][j]
            ok &= lp
            print(f'(a,b) k={k} w={[str(x) for x in w]} bscale={bsc}: all ten blocks equal j1lib entry by entry: {good}; '
                  f'Lemma P2 identities: {lp} [{time.time() - t0:.1f}s]', flush=True)


# ---------------- (c): two levels, reduced forms against projections
def projections(w, levels, bsc):
    """levels: list of level labels (0/1) of the copies 2..k-1.  Returns the projected invariant and standard forms of
    Phi (bsc = 1 only) and of the gauge-fixed block G."""
    k = len(w)
    st, M, Mt = jblock(w, bsc)
    n = len(st)
    R, Rp, P, Q = index_sets(st, k)
    orbit = lambda c: 't' if c == 0 else ('T' if c == 1 else levels[c - 2])
    out = {}
    # Phi = C - B^T diag(1/N) B (unrescaled), then rescaled
    if bsc == 1:
        Nr = {r: (1 if r in (0, 1) else 2) for r in range(k)}
        oth = R + Rp + Q
        Phi = [[M[i][j] - sum(M[l][i] * M[l][j] / Nr[st[l][0][1]] for l in P) for j in oth] for i in oth]
        mu = {i: w[copy_of(st[i])] for i in range(n)}
        Phit = [[Phi[x][y] * mu[oth[x]] * mu[oth[y]] for y in range(len(oth))] for x in range(len(oth))]
        out['Phi'] = (oth, Phit)
    keep = [i for i in range(n) if not (i in P and st[i][0][1] == 0)]
    G = [[Mt[i][j] for j in keep] for i in keep]
    out['G'] = (keep, G)
    res = {}
    for key, (idx, Mx) in out.items():
        # test vectors: for each state type (in the order of the algebra) and each orbit of its index, the indicator
        types = [R, Rp, P, Q]
        vecs = []
        for tp in types:
            orbs = []
            for i in idx:
                if i in tp:
                    o = orbit(copy_of(st[i]))
                    if o not in orbs:
                        orbs.append(o)
            order = [o for o in ('t', 'T', 0, 1) if o in orbs]
            for o in order:
                vecs.append([1 if (idx[x] in tp and orbit(copy_of(st[idx[x]])) == o) else 0 for x in range(len(idx))])
        F = [[sum(v1[x] * Mx[x][y] * v2[y] for x in range(len(idx)) for y in range(len(idx)) if v1[x] and v2[y])
              for v2 in vecs] for v1 in vecs]
        std = []
        for lv in (0, 1):
            hs = [c for c in range(2, k) if levels[c - 2] == lv]
            if len(hs) < 2:
                std.append(None)
                continue
            h1, h2 = hs[0], hs[1]
            fv = []
            for tp in types:
                if any(i in tp for i in idx):
                    fv.append([(1 if copy_of(st[idx[x]]) == h1 else (-1 if copy_of(st[idx[x]]) == h2 else 0))
                               if idx[x] in tp else 0 for x in range(len(idx))])
            S = [[sum(v1[x] * Mx[x][y] * v2[y] for x in range(len(idx)) for y in range(len(idx)) if v1[x] and v2[y]) / 2
                  for v2 in fv] for v1 in fv]
            std.append(S)
        res[key] = (F, std)
    return res


size_pairs = [(1, 1), (2, 1), (1, 2), (2, 2), (3, 1), (3, 2), (2, 3), (4, 1), (3, 3), (2, 0), (3, 0), (4, 0)]
for (u0, u1) in size_pairs:
    k = 2 + u0 + u1
    a = Fr(random.randint(1, 9), random.randint(1, 9))
    b = Fr(random.randint(1, 9), random.randint(1, 9))
    wq = Fr(random.randint(1, 9), random.randint(1, 9))
    levels = [0] * u0 + [1] * u1
    random.shuffle(levels)
    w = [a, b] + [Fr(1) if lv == 0 else wq for lv in levels]
    s = [x * x for x in w]
    for bsc in (Fr(1), Fr(1001, 1000)):
        proj = projections(w, levels, bsc)
        qs = [Fr(1), wq * wq] if u1 > 0 else [Fr(1)]
        us = [u0, u1] if u1 > 0 else [u0]
        alg = t2lev.point(s[0], s[1], qs, us, bsc)
        good = True
        for key in proj:
            Finv, Fstd = t2lev.forms(alg, key)
            Pinv, Pstd = proj[key]
            good &= Finv == Pinv
            for lv in range(len(us)):
                if Pstd[lv] is not None:
                    good &= Fstd[lv] == Pstd[lv]
        ok &= good
        print(f'(c) (u0,u1)=({u0},{u1}) a={s[0]} b={s[1]} q={wq * wq} bscale={bsc}: forms {sorted(proj)} equal the '
              f'projections: {good} [{time.time() - t0:.1f}s]', flush=True)
print('PASS' if ok else 'FAIL')
