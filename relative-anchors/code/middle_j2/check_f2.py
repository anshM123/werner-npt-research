"""Exact checks of the family-2 closed forms and reductions (f2sym.py) against the Lemma B builder j1lib.py (j = 2).

(a) Closed forms, entry by entry: at two-level rational weight vectors (concrete copies, every equality pattern), every
    entry of the blocks R-R, P-P, R-P, P-Q equals the rescaled j1lib.py block, and the blocks R-Q and Q-Q equal the
    formulas -E~ + 1/2 Ahat2^-1 Ncheck2 D^2 and 2 W4^-1 built as explicit matrices on all pairs; B-part scaled by 1
    and by 1001/1000.
(b) Reduced forms: at two-level points, for every part, the reduced blocks (RR, RP, RQ, PP, PQ, QQ) and the Gram
    matrices of f2sym.reduce_part2 equal the projections of the rescaled j1lib.py block on the explicit test vectors
    (named copies placed on concrete copies), the kernel vectors are annihilated, and the dimension counts add up to
    the block size.
Usage: python check_f2.py"""
import itertools, random, time
from fractions import Fraction as Fr
from j1lib import J1
from r2lib import states_T1
from f2sym import Model2, reduce_part2, part_spec2, keep_for2, assemble
from l1sym import minv, mm

t0 = time.time()
random.seed(23)
ok = True


def jblock(w, bscale):
    k = len(w)
    st = states_T1(k, 0)
    M = J1(w, bscale=bscale).matrix(st)
    mu = []
    for (a, y) in st:
        m = Fr(1)
        if a is not None and a[0] == 1:
            m *= w[a[1]]
        for c in range(k):
            if y[c] == 1:
                m *= w[c]
        mu.append(m)
    Mt = [[M[i][j] * mu[i] * mu[j] for j in range(len(st))] for i in range(len(st))]
    return st, Mt


def concrete_labels(levels, named):
    """levels: level (0/1) of the copies 1..k-1 (copy 0 is t).  named: the named labels of a part; each is placed on
    the first free copy of its level.  Returns the list of labels of the copies."""
    labs = [(2, 't', 't')]
    used = set()
    for i, lv in enumerate(levels):
        labs.append((lv, None, ('c', i)))
    for nl in named:
        if nl[0] == 2:
            continue
        j = next(c for c in range(1, len(labs)) if labs[c][0] == nl[0] and c not in used)
        used.add(j)
        labs[j] = nl
    return labs


def kind(state):
    a, y = state
    if a == (2, 0):
        return 'R'
    if a is None:
        return 'Q'
    return 'P'


def idx_of(state):
    a, y = state
    ones = [c for c in range(len(y)) if y[c] == 1]
    if kind(state) == 'P':
        return (a[1], ones[0])
    return tuple(ones)


def model_at(w, levels, bsc):
    s = [x * x for x in w]
    q = [None, None]
    for i, lv in enumerate(levels):
        q[lv] = s[i + 1]
    u = [levels.count(0), levels.count(1)]
    Q = [q[0] if q[0] is not None else Fr(7, 3), q[1] if q[1] is not None else Fr(5, 2)]
    return Model2(Fr(1), s[0], Q, u, bsc), u


# ---------------- (a) entries
for (u0, u1) in [(2, 1), (1, 2), (2, 2), (3, 1), (2, 3), (3, 0)]:
    levels = [0] * u0 + [1] * u1
    random.shuffle(levels)
    k = 1 + u0 + u1
    wt = Fr(random.randint(1, 9), random.randint(1, 9))
    wl = [Fr(random.randint(1, 9), random.randint(1, 9)), Fr(random.randint(1, 9), random.randint(1, 9))]
    w = [wt] + [wl[lv] for lv in levels]
    for bsc in (Fr(1), Fr(1001, 1000)):
        st, Mt = jblock(w, bsc)
        M, u = model_at(w, levels, bsc)
        labs = concrete_labels(levels, [])
        n = len(st)
        good = True
        Rst = [i for i in range(n) if kind(st[i]) == 'R']
        Pst = [i for i in range(n) if kind(st[i]) == 'P']
        Qst = [i for i in range(n) if kind(st[i]) == 'Q']
        for i in range(n):
            for j in range(n):
                ki, kj = kind(st[i]), kind(st[j])
                ii, jj = idx_of(st[i]), idx_of(st[j])
                L1 = [labs[c] for c in ii]
                L2 = [labs[c] for c in jj]
                if ki == 'R' and kj == 'R':
                    v = M.RR(tuple(L1), tuple(L2))
                elif ki == 'P' and kj == 'P':
                    v = M.PP(L1[0], L1[1], L2[0], L2[1])
                elif ki == 'R' and kj == 'P':
                    v = M.RP(L1[0], L1[1], L2[0], L2[1])
                elif ki == 'P' and kj == 'Q':
                    v = M.PQ(L1[0], L1[1], L2[0], L2[1])
                else:
                    continue
                if v != Mt[i][j]:
                    good = False
        # R-Q and Q-Q through explicit pair matrices
        allpairs = list(itertools.combinations(range(k), 2))
        Upairs = [P for P in allpairs if 0 not in P]
        lab = lambda P: (labs[P[0]], labs[P[1]])
        A2 = [[M.A2hat(*lab(P), *lab(P2)) for P2 in allpairs] for P in allpairs]
        N2 = [[M.N2c(*lab(P), *lab(P2)) for P2 in Upairs] for P in allpairs]
        sP = lambda P: M.s(labs[P[0]]) * M.s(labs[P[1]])
        A2i = minv(A2, Fr(1))
        T = mm(A2i, N2, Fr(0))
        RQ = {(P, P2): -(sP(P) if P == P2 else 0) + bsc * T[a][b] * sP(P2) / 2
              for a, P in enumerate(allpairs) for b, P2 in enumerate(Upairs)}
        W = [[M.W4(*lab(P), *lab(P2)) for P2 in Upairs] for P in Upairs]
        Wi = minv(W, Fr(1))
        QQ = {(P, P2): 2 * (sP(P) if P == P2 else 0) + bsc * (2 * Wi[a][b] - 2 * (sP(P) if P == P2 else 0))
              for a, P in enumerate(Upairs) for b, P2 in enumerate(Upairs)}
        for i in Rst:
            for j in Qst:
                good &= RQ[(idx_of(st[i]), idx_of(st[j]))] == Mt[i][j]
        for i in Qst:
            for j in Qst:
                good &= QQ[(idx_of(st[i]), idx_of(st[j]))] == Mt[i][j]
        ok &= good
        print(f'(a) (u0,u1)=({u0},{u1}) w={[str(x) for x in w]} bscale={bsc}: all six blocks equal j1lib: {good} '
              f'[{time.time() - t0:.1f}s]', flush=True)


# ---------------- (b) reduced forms
def project(st, Mt, labs, part, u0, u1):
    named, R, P, Qv, ker0, gauge0, names = part_spec2(part)
    kR, kP, kQ, gauge, ker = keep_for2(part, min(u0, 3), min(u1, 3))
    vec = {'R': [R[i] for i in kR], 'P': [P[i] for i in kP], 'Q': [Qv[i] for i in kQ]}
    n = len(st)
    cols = {}
    for key in 'RPQ':
        cols[key] = []
        for f in vec[key]:
            v = [0] * n
            for i in range(n):
                if kind(st[i]) == key:
                    ii = idx_of(st[i])
                    v[i] = f(labs[ii[0]], labs[ii[1]])
            cols[key].append(v)
    nz = lambda v: [i for i in range(n) if v[i]]
    proj = lambda v1, v2: sum(v1[i] * Mt[i][j] * v2[j] for i in nz(v1) for j in nz(v2))
    out = {}
    for a, b in [('R', 'R'), ('R', 'P'), ('R', 'Q'), ('P', 'P'), ('P', 'Q'), ('Q', 'Q')]:
        out[a + b] = [[proj(v1, v2) for v2 in cols[b]] for v1 in cols[a]]
    for key in 'RPQ':
        out['G' + key] = [[sum(x * y for x, y in zip(v1, v2)) for v2 in cols[key]] for v1 in cols[key]]
    return out, cols


def dims(part, u0, u1):
    """dimension of the irrep (u0, u1 >= the needed sizes)."""
    std = lambda u: u - 1
    return {'inv': 1, 'std': std(u0), 'two': u0 * (u0 - 3) // 2, 'anti': (u0 - 1) * (u0 - 2) // 2,
            'mix': std(u0) * std(u1)}[part]


def parts_at(u0, u1):
    out = [('inv', False)]
    for sw in (False, True):
        a, b = (u1, u0) if sw else (u0, u1)
        if a >= 2:
            out.append(('std', sw))
        if a >= 4:
            out.append(('two', sw))
        if a >= 3:
            out.append(('anti', sw))
    if u0 >= 2 and u1 >= 2:
        out.append(('mix', False))
    return out


for (u0, u1) in [(2, 1), (1, 2), (2, 2), (3, 1), (1, 3), (3, 2), (2, 3), (4, 1), (4, 2), (2, 4), (3, 3),
                 (4, 0), (5, 0)]:
    k = 1 + u0 + u1
    wt = Fr(random.randint(1, 9), random.randint(1, 9))
    wl = [Fr(random.randint(1, 9), random.randint(1, 9)), Fr(random.randint(1, 9), random.randint(1, 9))]
    levels = [0] * u0 + [1] * u1
    random.shuffle(levels)
    w = [wt] + [wl[lv] for lv in levels]
    for bsc in (Fr(1), Fr(1001, 1000)):
        st, Mt = jblock(w, bsc)
        good = True
        cnt = {'R': 0, 'P': 0, 'Q': 0}
        for part, sw in parts_at(u0, u1):
            # sw: the irrep on level 1: exchange the level labels
            lv2 = [1 - x for x in levels] if sw else levels
            uu = (u1, u0) if sw else (u0, u1)
            wl2 = (wl[1], wl[0]) if sw else (wl[0], wl[1])
            M = Model2(Fr(1), wt * wt, [wl2[0] ** 2, wl2[1] ** 2], list(uu), bsc)
            red = reduce_part2(M, part, min(uu[0], 3), min(uu[1], 3))
            named = part_spec2(part)[0]
            labs = concrete_labels(lv2, named)
            pr, cols = project(st, Mt, labs, part, uu[0], uu[1])
            for key in ['RR', 'RP', 'RQ', 'PP', 'PQ', 'QQ', 'GR', 'GP', 'GQ']:
                if key in red:
                    good &= red[key] == pr[key]
                else:
                    good &= all(len(r) == 0 for r in pr[key]) or len(pr[key]) == 0
            if bsc == 1:
                F, iR, iP, iQ = assemble(red)
                for kv in red['ker']:
                    vec = [0] * len(F)
                    for j, c in enumerate(kv):
                        vec[iP[j]] = c
                    good &= all(sum(F[i][j] * vec[j] for j in range(len(F))) == 0 for i in range(len(F)))
            d = dims(part, *uu)
            for key in 'RPQ':
                cnt[key] += d * len(red['names']['RPQ'.index(key)])
        expect = {'R': k * (k - 1) // 2, 'P': k * (k - 1), 'Q': (k - 1) * (k - 2) // 2}
        good &= cnt == expect
        ok &= good
        print(f'(b) (u0,u1)=({u0},{u1}) bscale={bsc}: parts {[p + ("[1]" if s else "") for p, s in parts_at(u0, u1)]}: '
              f'reduced blocks = projections, kernels, dimensions {cnt}: {good} [{time.time() - t0:.1f}s]', flush=True)
print('PASS' if ok else 'FAIL')
