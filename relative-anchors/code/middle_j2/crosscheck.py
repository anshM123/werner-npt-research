"""Second, independent checks of the (R_inf)_2 results (families 1 and 2) against the builder from the definitions
(ylib.py), written separately from t2lev.py and f2sym.py (only their certified forms are evaluated here).

A. Family 1, exact: at two-level rational points (k = 7, 8; sizes and placements other than those of check_t2.py), the block is
   built by ylib (rescaled frame, units 2^{k-2}); P_t is deleted, the exact Schur complement onto (R, R', Q) is taken
   (this is Phi~ of Lemma P2, computed without the diag(1/N) formula), projected on the orbit indicators and on
   e_h - e_h' of each level; these must equal the symbolic forms that the certificates use (t2lev.symbolic in the
   certified regime), evaluated at the point.
B. Family 2, spectral: the eigenvalues of the full rescaled ylib block must equal the union of the generalized
   eigenvalues of the reduced forms (F, Gram) of all parts, each with the dimension of its irrep as multiplicity.
   The reduced forms are the symbolic ones of the certified regimes (parts on level 0; a part on level 1 through the
   level exchange (S1): sizes exchanged, t-weight level flipped, eigenvalues times q^2).
C. End to end, exact: kernel vectors (column-group indicators in the rescaled frame) are in the kernel, the nullity is
   their number, and the block with one P-state per column group deleted is positive definite (exact LDL).
Usage: python crosscheck.py [A|B|C|all]"""
import sys, time, random
from fractions import Fraction as Fr

WHICH = sys.argv[1] if len(sys.argv) > 1 else 'all'
import flint
import numpy as np
import scipy.linalg as sla
from ylib import Builder, point_ring
from j1lib import to_ylib_state
from r2lib import states_T2, states_T1
import t2lev
from flintfield import Field
from f2sym import Model2, reduce_part2, assemble

t0 = time.time()
ok = True
random.seed(20261105)


def fq(x):
    return flint.fmpq(x.numerator, x.denominator)


def ylib_block(s, states):
    k = len(s)
    Bd = Builder(point_ring([fq(x) for x in s]))
    ys = [to_ylib_state(x, k) for x in states]
    n = len(states)
    M = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            v = Bd.Y(ys[i], ys[j])
            M[i][j] = M[j][i] = Fr(int(v.p), int(v.q)) * 2 ** (k - 2)
    return M


def ev(f, names, vals):
    """evaluate a flintfield element at the point vals (dict name -> Fraction)."""
    args = [fq(vals[nm]) for nm in names]
    if len(names) == 0:
        n = f.n.leading_coefficient() if not f.n.is_zero() else flint.fmpq(0)
        d = f.d.leading_coefficient()
    else:
        n, d = f.n(*args), f.d(*args)
    n, d = Fr(int(n.p), int(n.q)), Fr(int(d.p), int(d.q))
    return n / d


def solve(A, B):
    """A^{-1} B exactly (A square invertible, B a list of columns as a matrix)."""
    n = len(A)
    m = len(B[0])
    M = [list(A[i]) + list(B[i]) for i in range(n)]
    for c in range(n):
        r = next(i for i in range(c, n) if M[i][c] != 0)
        M[c], M[r] = M[r], M[c]
        pv = M[c][c]
        M[c] = [v / pv for v in M[c]]
        for i in range(n):
            if i != c and M[i][c] != 0:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    return [row[n:] for row in M]


def ldl_pd(A):
    A = [list(r) for r in A]
    n = len(A)
    for c in range(n):
        p = A[c][c]
        if p <= 0:
            return False
        for i in range(c + 1, n):
            if A[i][c] != 0:
                f = A[i][c] / p
                for j in range(c + 1, n):
                    A[i][j] -= f * A[c][j]
    return True


def rank(A):
    return flint.fmpq_mat([[fq(x) for x in r] for r in A]).rank()


def end_to_end(M, states, label):
    """C: kernel vectors, nullity, gauge-fixed PD (rescaled frame: kernel vectors are column-group indicators)."""
    global ok
    groups = {}
    for i, (a, y) in enumerate(states):
        if a is not None and a[0] == 1:
            groups.setdefault(y, []).append(i)
    n = len(states)
    inker = True
    for y, idx in groups.items():
        v = [0] * n
        for i in idx:
            v[i] = 1
        inker &= all(sum(M[i][j] for j in idx) == 0 for i in range(n))
    null = n - rank(M)
    drop = set(idx[-1] for idx in groups.values())
    keep = [i for i in range(n) if i not in drop]
    pd = ldl_pd([[M[i][j] for j in keep] for i in keep])
    good = inker and null == len(groups) and pd
    ok &= good
    print(f'  C {label}: size {n}, column groups {len(groups)} in kernel: {inker}, nullity {null}, gauge-fixed PD: {pd}'
          f' -> {good} [{time.time() - t0:.1f}s]', flush=True)


# ------------------------------------------------------------------ A: family 1
def fam1(levels, tw, sv, spec):
    """levels: level of each U-copy (copies 2..k-1); tw = (level of t, level of tau) or None = free weight;
    sv = dict with 'a', 'b' (if free) and 'q'; spec = list of regime specs per level ('2+' or fixed)."""
    global ok
    l = len(spec)
    k = 2 + len(levels)
    qv = [Fr(1), sv.get('q', Fr(1))]
    s = [sv['a'] if tw[0] is None else qv[tw[0]], sv['b'] if tw[1] is None else qv[tw[1]]] + [qv[x] for x in levels]
    st = states_T2(k, 0, 1)
    M = ylib_block(s, st)
    label = f'family 1, k={k}, U-levels {levels}, T-weights {tw}, s={[str(x) for x in s]}'
    if WHICH in ('C', 'all'):
        end_to_end(M, st, label)
    if WHICH not in ('A', 'all'):
        return
    iP = [i for i, (a, y) in enumerate(st) if a is not None and a[0] == 1 and a[1] != 0]
    iC = [i for i, (a, y) in enumerate(st) if not (a is not None and a[0] == 1)]
    Pi = [[M[i][j] for j in iP] for i in iP]
    Bm = [[M[i][j] for j in iC] for i in iP]
    pdP = ldl_pd(Pi)
    X = solve(Pi, Bm)
    Phi = [[M[iC[a]][iC[b]] - sum(Bm[r][a] * X[r][b] for r in range(len(iP))) for b in range(len(iC))]
           for a in range(len(iC))]
    cst = [st[i] for i in iC]

    def copy_of(state):
        a, y = state
        return [c for c in range(k) if y[c] == 1][0]

    def orb(c):
        return 't' if c == 0 else ('T' if c == 1 else levels[c - 2])
    kinds = [lambda x: x[0] == (2, 0), lambda x: x[0] == (2, 1), lambda x: x[0] is None]
    spaces = [('t',) + tuple(range(l)), ('T',) + tuple(range(l)), tuple(range(l))]
    vecs = []
    for kd, sp in zip(kinds, spaces):
        for o in sp:
            vecs.append([1 if (kd(x) and orb(copy_of(x)) == o) else 0 for x in cst])
    mine_inv = [[sum(v1[i] * Phi[i][j] * v2[j] for i in range(len(cst)) if v1[i] for j in range(len(cst)) if v2[j])
                 for v2 in vecs] for v1 in vecs]
    fixed = [None if sp.endswith('+') else int(sp) for sp in spec]
    offs = [int(sp[:-1]) if sp.endswith('+') else 0 for sp in spec]
    alg, K, names = t2lev.symbolic(l, fixed=fixed, offsets=offs, tw=tw)
    Finv, Fstd = t2lev.forms(alg, 'Phi')
    u = [levels.count(i) for i in range(l)]
    vals = {'a': sv.get('a'), 'b': sv.get('b'), 'q1': qv[1]}
    for i in range(l):
        if fixed[i] is None:
            vals[f'x{i}'] = Fr(u[i] - offs[i])
        else:
            assert fixed[i] == u[i]
    theirs_inv = [[ev(f, names, vals) for f in row] for row in Finv]
    good = theirs_inv == mine_inv
    for lv in range(l):
        hs = [c for c in range(2, k) if levels[c - 2] == lv]
        if len(hs) < 2:
            continue
        fv = [[(1 if copy_of(x) == hs[0] else (-1 if copy_of(x) == hs[1] else 0)) if kd(x) else 0 for x in cst]
              for kd in kinds]
        mine_std = [[sum(v1[i] * Phi[i][j] * v2[j] for i in range(len(cst)) if v1[i] for j in range(len(cst)) if v2[j])
                     / 2 for v2 in fv] for v1 in fv]
        theirs_std = [[ev(f, names, vals) for f in row] for row in Fstd[lv]]
        good &= theirs_std == mine_std
    ok &= good and pdP
    print(f'  A {label}, regime {spec}: Pi_g PD {pdP}; my Schur complement projected = certified symbolic forms: {good}'
          f' [{time.time() - t0:.1f}s]', flush=True)


# ------------------------------------------------------------------ B: family 2
REG = {('inv', 0): [('1', '3+'), ('3+', '1'), ('2', '2'), ('2', '3+'), ('3+', '2')] +
       [(str(m), '3+') for m in range(3, 13)] + [('3+', str(m)) for m in range(3, 13)],
       ('mix', 0): [('2', '2'), ('2', '3+'), ('3+', '2'), ('3+', '3+')]}
for lt in (0, 1):
    REG[('std', lt)] = [('2', '2'), ('2', '3+'), ('3+', '1'), ('3+', '2'), ('3+', '3+')]
    REG[('two', lt)] = [('4+', '1'), ('4+', '2'), ('4+', '3+')]
    REG[('anti', lt)] = [('3+', '1'), ('3+', '2'), ('3+', '3+')]
ONE = {'inv': '3+', 'std': '3+', 'two': '4+', 'anti': '3+'}
DIM = {'inv': lambda a, b: 1, 'std': lambda a, b: a - 1, 'two': lambda a, b: a * (a - 3) // 2,
       'anti': lambda a, b: (a - 1) * (a - 2) // 2, 'mix': lambda a, b: (a - 1) * (b - 1)}
NEED = {'inv': lambda a, b: True, 'std': lambda a, b: a >= 2, 'two': lambda a, b: a >= 4,
        'anti': lambda a, b: a >= 3, 'mix': lambda a, b: a >= 2 and b >= 2}


def inreg(sp, u):
    return u >= int(sp[:-1]) if sp.endswith('+') else u == int(sp)


def part_eigs(part, tw, v0, v1, a, q):
    """generalized eigenvalues of the certified reduced form of part (level 0) at sizes (v0, v1), t-weight level tw
    (or tw = 'free' with squared weight a; one level), level weights (1, q)."""
    if v1 == 0:
        spec = (ONE[part], '0')
    else:
        cand = [r for r in REG[(part, tw)] if inreg(r[0], v0) and inreg(r[1], v1)]
        spec = cand[0]
    f = [None if sp.endswith('+') else int(sp) for sp in spec]
    o = [int(sp[:-1]) if sp.endswith('+') else None for sp in spec]
    names = (['a'] if tw == 'free' else []) + (['q'] if (f[1] is None or f[1] > 0 or tw == 1) else []) + \
        (['x0'] if f[0] is None else []) + (['x1'] if f[1] is None else [])
    K = Field(names)
    g = dict(zip(names, K.gens))
    one = K(1)
    qq = g.get('q', K(5) / K(2))
    U = [g['x0'] + o[0] if f[0] is None else K(f[0]), g['x1'] + o[1] if f[1] is None else K(f[1])]
    aa = g['a'] if tw == 'free' else (one if tw == 0 else qq)
    Mod = Model2(one, aa, [one, qq], U)
    red = reduce_part2(Mod, part, f[0] if f[0] is not None else o[0], f[1] if f[1] is not None else o[1])
    F, iR, iP, iQ = assemble(red)
    n = len(F)
    G = [[K(0)] * n for _ in range(n)]
    for key, idx in (('GR', iR), ('GP', iP), ('GQ', iQ)):
        if idx:
            for i in range(len(idx)):
                for j in range(len(idx)):
                    G[idx[i]][idx[j]] = red[key][i][j]
    vals = {'a': a, 'q': q, 'x0': Fr(v0 - (o[0] or 0)), 'x1': Fr(v1 - (o[1] or 0))}
    Fn = np.array([[float(ev(x, names, vals)) for x in r] for r in F])
    Gn = np.array([[float(ev(x, names, vals)) for x in r] for r in G])
    return list(sla.eigh(Fn, Gn, eigvals_only=True)), spec


def fam2(levels, tw, a, q):
    """levels of the U-copies (copies 1..k-1); t = copy 0 with the weight of level tw, or tw = 'free' (weight a)."""
    global ok
    k = 1 + len(levels)
    qv = [Fr(1), q]
    s = [a if tw == 'free' else qv[tw]] + [qv[x] for x in levels]
    st = states_T1(k, 0)
    M = ylib_block(s, st)
    label = f'family 2, k={k}, U-levels {levels}, t-weight {tw}, s={[str(x) for x in s]}'
    if WHICH in ('C', 'all'):
        end_to_end(M, st, label)
    if WHICH not in ('B', 'all'):
        return
    full = sorted(np.linalg.eigvalsh(np.array([[float(x) for x in r] for r in M])))
    u0, u1 = levels.count(0), levels.count(1)
    union = []
    used = []
    for part in ('inv', 'mix', 'std', 'two', 'anti'):
        for pl in ((None,) if part in ('inv', 'mix') else (0, 1)):
            sizes = (u0, u1) if pl in (None, 0) else (u1, u0)
            if not NEED[part](*sizes):
                continue
            if part == 'mix' and (u0 == 0 or u1 == 0):
                continue
            mult = DIM[part](*sizes)
            if mult == 0:
                continue
            direct = pl in (None, 0) and (u1 == 0 or any(inreg(r[0], u0) and inreg(r[1], u1)
                                                          for r in REG.get((part, tw), [])))
            if direct:
                e, spec = part_eigs(part, tw, u0, u1, a, q)
                fac = 1
            else:            # (S1): exchange the levels, divide all squared weights by q (entries scale by 1/q^2)
                assert tw != 'free'
                e, spec = part_eigs(part, 1 - tw, u1, u0, None, 1 / q)
                spec = ('S1',) + tuple(spec)
                fac = float(q) ** 2
            union += [x * fac for x in e] * mult
            used.append(f'{part}{"" if pl is None else pl}{spec}x{mult}')
    union = sorted(union)
    same_n = len(union) == len(full)
    err = max(abs(x - y) for x, y in zip(union, full)) if same_n else float('inf')
    scale = max(abs(x) for x in full)
    good = same_n and err <= 1e-9 * scale
    ok &= good
    print(f'  B {label}: {len(full)} eigenvalues; parts {used}; union of part spectra = full spectrum: {good} '
          f'(max deviation {err:.2e}, scale {scale:.2e}) [{time.time() - t0:.1f}s]', flush=True)


if WHICH in ('A', 'C', 'all'):
    print('family 1', flush=True)
    fam1([0, 0, 0, 0, 0], (None, None), {'a': Fr(7, 3), 'b': Fr(2, 9)}, ['2+'])                    # F1a, k = 7
    fam1([0, 1, 1, 0, 1, 1], (0, 0), {'q': Fr(1, 40)}, ['2+', '2+'])                               # k = 8
    fam1([1, 0, 1, 1, 1, 1], (0, 0), {'q': Fr(40)}, ['1', '2+'])                                   # k = 8
    fam1([0, 0, 1, 0, 0, 0], (0, 0), {'q': Fr(5, 3)}, ['2+', '1'])                                 # k = 8
    fam1([1, 0, 0, 1, 1, 0], (0, 1), {'q': Fr(1, 900)}, ['3', '2+'])                               # k = 8, m = 3
    fam1([0, 1, 1, 1, 1, 0], (0, 1), {'q': Fr(30)}, ['2', '2+'])                                   # k = 8, m = 2
    fam1([1, 1, 0, 1, 1], (0, 1), {'q': Fr(2, 7)}, ['1', '2+'])                                    # k = 7
if WHICH in ('B', 'C', 'all'):
    print('family 2', flush=True)
    fam2([0, 0, 0, 0, 0, 0], 'free', Fr(1, 13), Fr(1))                                             # F2a, k = 7
    fam2([0, 1, 1, 0, 1, 0, 1], 0, None, Fr(1, 30))                                                # (3, 4), k = 8
    fam2([1, 1, 0, 1, 1, 1, 0], 1, None, Fr(17, 5))                                                # (2, 5), k = 8
    fam2([0, 0, 0, 0, 1, 0, 0], 1, None, Fr(1, 900))                                               # (6, 1), k = 8
    fam2([0, 1, 0, 1, 0, 1], 0, None, Fr(25))                                                      # (3, 3), k = 7
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
