"""The reduced forms of families 1 and 2 against the builder from the definitions (ylib.py) directly, at two-level
rational weight vectors: the rescaled block is built by ylib.Builder (M, V_p, D_p, R_p and the partial transpose; not
Lemma B), multiplied by 2^{k-2}, and projected on the explicit test vectors.  The projections must equal, exactly,
the reduced forms of t2lev.py (invariant and standard forms of the gauge-fixed family-1 block, P_t deleted) and of
f2sym.py (every part, every reduced block of family 2).  Usage: python check_ylib_reduced.py"""
import time
from fractions import Fraction as Fr
import flint
from ylib import Builder, point_ring
from j1lib import to_ylib_state
from r2lib import states_T2, states_T1
from f2sym import Model2, reduce_part2, part_spec2, keep_for2
import t2lev

t0 = time.time()
ok = True


def ylib_block(w, states):
    k = len(w)
    Bd = Builder(point_ring([flint.fmpq(x.numerator ** 2, x.denominator ** 2) for x in w]))
    ys = [to_ylib_state(x, k) for x in states]
    n = len(states)
    M = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            v = Bd.Y(ys[i], ys[j])
            M[i][j] = M[j][i] = Fr(int(v.p), int(v.q)) * 2 ** (k - 2)
    return M


def copy_of(state):
    a, y = state
    if a is not None and a[0] == 1:
        return a[1]
    return [c for c in range(len(y)) if y[c] == 1][0]


def quad(v1, M, v2):
    n = len(M)
    return sum(v1[i] * M[i][j] * v2[j] for i in range(n) if v1[i] for j in range(n) if v2[j])


# ---------------- family 1 (copies 0 = t, 1 = tau, then U with the given levels)
for (levels, w) in [([0, 0, 1], [Fr(1, 2), Fr(3), Fr(1), Fr(1), Fr(2)]),
                    ([0, 1, 1], [Fr(2), Fr(1), Fr(1), Fr(1, 3), Fr(1, 3)]),
                    ([0, 1, 0, 1], [Fr(1), Fr(1, 2), Fr(1), Fr(3), Fr(1), Fr(3)])]:
    k = len(w)
    st = states_T2(k, 0, 1)
    M = ylib_block(w, st)
    keep = [i for i in range(len(st)) if not (st[i][0] is not None and st[i][0][0] == 1 and st[i][0][1] == 0)]
    G = [[M[i][j] for j in keep] for i in keep]
    sts = [st[i] for i in keep]
    orbit = lambda c: 't' if c == 0 else ('T' if c == 1 else levels[c - 2])
    types = [lambda s: s[0] == (2, 0), lambda s: s[0] == (2, 1),
             lambda s: s[0] is not None and s[0][0] == 1, lambda s: s[0] is None]
    vecs = []
    for tp in types:
        orbs = []
        for s in sts:
            if tp(s) and orbit(copy_of(s)) not in orbs:
                orbs.append(orbit(copy_of(s)))
        for o in [o for o in ('t', 'T', 0, 1) if o in orbs]:
            vecs.append([1 if (tp(s) and orbit(copy_of(s)) == o) else 0 for s in sts])
    P_inv = [[quad(v1, G, v2) for v2 in vecs] for v1 in vecs]
    u0, u1 = levels.count(0), levels.count(1)
    s = [x * x for x in w]
    qv = [s[2 + levels.index(0)], s[2 + levels.index(1)]]
    alg = t2lev.point(s[0], s[1], qv, [u0, u1])
    Finv, Fstd = t2lev.forms(alg, 'G')
    good = Finv == P_inv
    for lv in (0, 1):
        hs = [c for c in range(2, k) if levels[c - 2] == lv]
        if len(hs) < 2:
            continue
        fv = []
        for tp in types:
            if any(tp(x) for x in sts):
                fv.append([(1 if copy_of(x) == hs[0] else (-1 if copy_of(x) == hs[1] else 0)) if tp(x) else 0
                           for x in sts])
        S = [[quad(v1, G, v2) / 2 for v2 in fv] for v1 in fv]
        good &= S == Fstd[lv]
    ok &= good
    print(f'family 1, levels {levels}, w = {[str(x) for x in w]}: gauge-fixed reduced forms = projections of the ylib '
          f'block: {good} [{time.time() - t0:.1f}s]', flush=True)


# ---------------- family 2 (copy 0 = t, then U with the given levels)
def labels_for(levels, named):
    labs = [(2, 't', 't')] + [(lv, None, ('c', i)) for i, lv in enumerate(levels)]
    used = set()
    for nl in named:
        if nl[0] == 2:
            continue
        j = next(c for c in range(1, len(labs)) if labs[c][0] == nl[0] and c not in used)
        used.add(j)
        labs[j] = nl
    return labs


def kind(state):
    a, y = state
    return 'R' if a == (2, 0) else ('Q' if a is None else 'P')


def idx(state):
    a, y = state
    ones = [c for c in range(len(y)) if y[c] == 1]
    return (a[1], ones[0]) if kind(state) == 'P' else tuple(ones)


for (levels, w) in [([0, 0, 1, 1], [Fr(1, 2), Fr(1), Fr(1), Fr(3), Fr(3)]),
                    ([0, 0, 0, 1], [Fr(2), Fr(1), Fr(1), Fr(1), Fr(1, 3)]),
                    ([0, 1, 0, 1, 1], [Fr(1), Fr(1), Fr(1, 2), Fr(1), Fr(1, 2), Fr(1, 2)])]:
    k = len(w)
    st = states_T1(k, 0)
    M = ylib_block(w, st)
    u0, u1 = levels.count(0), levels.count(1)
    wl = [w[1 + levels.index(0)], w[1 + levels.index(1)]]
    good = True
    parts = [('inv', False)] + [(p, sw) for sw in (False, True) for p, need in
                                (('std', 2), ('two', 4), ('anti', 3)) if ((u1 if sw else u0) >= need)]
    if u0 >= 2 and u1 >= 2:
        parts.append(('mix', False))
    for part, sw in parts:
        lv2 = [1 - x for x in levels] if sw else levels
        uu = (u1, u0) if sw else (u0, u1)
        ww = (wl[1], wl[0]) if sw else (wl[0], wl[1])
        Mod = Model2(Fr(1), w[0] ** 2, [ww[0] ** 2, ww[1] ** 2], list(uu))
        red = reduce_part2(Mod, part, min(uu[0], 3), min(uu[1], 3))
        named, R, P, Qv, k0, g0, nm = part_spec2(part)
        kR, kP, kQ, gauge, ker = keep_for2(part, min(uu[0], 3), min(uu[1], 3))
        labs = labels_for(lv2, named)
        vec = {'R': [R[i] for i in kR], 'P': [P[i] for i in kP], 'Q': [Qv[i] for i in kQ]}
        cols = {}
        for key in 'RPQ':
            cols[key] = [[f(labs[idx(s)[0]], labs[idx(s)[1]]) if kind(s) == key else 0 for s in st] for f in vec[key]]
        for a, b in [('R', 'R'), ('R', 'P'), ('R', 'Q'), ('P', 'P'), ('P', 'Q'), ('Q', 'Q')]:
            if (a + b) in red:
                good &= red[a + b] == [[quad(v1, M, v2) for v2 in cols[b]] for v1 in cols[a]]
    ok &= good
    print(f'family 2, levels {levels}, w = {[str(x) for x in w]}: every reduced block of every part = projections of '
          f'the ylib block: {good} [{time.time() - t0:.1f}s]', flush=True)
print('PASS' if ok else 'FAIL')
