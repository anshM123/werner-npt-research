"""Explicit objects for the type I blocks with one marked copy, independent of the orbit calculus of tlib.py:
the states of a type I torus block of (R_inf)_j, the block built by fastj1.FastJ1 (Lemma B, exact) and rescaled by
mu (letters 1 in row and column), the explicit twisted vectors of tlib.NAMES as vectors over the states, the class sums
of the transpositions of U' and of T', and the known kernel vectors p (x) y of the block.

Copies: M = 0.  Case MT: T = {0, ..., tau - 1}, T' = {1, ..., tau - 1}, U' = {tau, ..., k - 1}.
Case MU: T = T' = {1, ..., tau}, U' = {tau + 1, ..., k - 1}.  Pairs: consecutive copies at the start of U'; R the
rest.  Anti mode: Ta, Tb the first two copies of T'."""
import itertools
from fractions import Fraction as Fr
from fastj1 import FastJ1
from j1lib import mu


def layout(k, tau, case, i):
    if case == 'MT':
        T = list(range(0, tau))
        Tp = list(range(1, tau))
        Up = list(range(tau, k))
    else:
        T = list(range(1, tau + 1))
        Tp = list(T)
        Up = list(range(tau + 1, k))
    pairs = [(Up[2 * l], Up[2 * l + 1]) for l in range(i)]
    R = Up[2 * i:]
    return T, Tp, Up, pairs, R


def block_states(k, j, T):
    T = sorted(T)
    tau = len(T)
    nonT = [c for c in range(k) if c not in T]
    st = []
    if j - tau >= 0:
        for r in range(k):
            for Y in itertools.combinations(nonT, j - tau):
                y = [0] * k
                for t in T:
                    y[t] = 2
                for c in Y:
                    y[c] = 1
                st.append(((1, r), tuple(y)))
    for t in T:
        free = sorted(nonT + [t])
        if 0 <= j - tau + 1 <= len(free):
            for Y in itertools.combinations(free, j - tau + 1):
                y = [0] * k
                for t2 in T:
                    if t2 != t:
                        y[t2] = 2
                for c in Y:
                    y[c] = 1
                st.append(((2, t), tuple(y)))
    if 0 <= j + 1 - tau <= len(nonT) and j + 1 <= k:
        for Y in itertools.combinations(nonT, j + 1 - tau):
            y = [0] * k
            for t in T:
                y[t] = 2
            for c in Y:
                y[c] = 1
            st.append((None, tuple(y)))
    return st


def block(k, j, T, w0, bscale=1):
    """rescaled block D (2^Z Y) D at weights (w0, 1, ..., 1), D = diag(mu)"""
    w = [Fr(w0)] + [Fr(1)] * (k - 1)
    st = block_states(k, j, T)
    M = FastJ1(w, bscale=bscale).matrix(st)
    mus = [mu(s, w) for s in st]
    n = len(st)
    return st, [[M[a][b] * mus[a] * mus[b] for b in range(n)] for a in range(n)]


def twisted(k, j, tau, case, mode, i, st):
    """dict name -> {state index: coefficient} for the explicit twisted vectors"""
    T, Tp, Up, pairs, R = layout(k, tau, case, i)
    pair_of = {}
    for l, (p, q) in enumerate(pairs):
        pair_of[p] = (l, q, 1)
        pair_of[q] = (l, p, -1)
    Ta, Tb = (Tp[0], Tp[1]) if len(Tp) >= 2 else (None, None)
    vecs = {}
    for idx, (a, y) in enumerate(st):
        skip = pair_of[a[1]][0] if (a is not None and a[1] in pair_of) else None
        sg = 1
        okp = True
        for l, (p, q) in enumerate(pairs):
            if l == skip:
                continue
            if y[p] + y[q] != 1:
                okp = False
                break
            sg *= 1 if y[p] == 1 else -1
        if not okp:
            continue
        name = None
        coef = sg
        if mode == 'anti':
            if a is None or a[1] not in (Ta, Tb):
                continue
            sgn = 1 if a[1] == Ta else -1
            coef *= sgn
            if case == 'MT':
                name = 'Ca' if a[0] == 1 else 'Da%d' % y[a[1]]
            else:
                name = ('Ca%d' % y[0]) if a[0] == 1 else 'Da%d%d' % (y[a[1]], y[0])
        else:
            if case == 'MT':
                if a is None:
                    name = 'G'
                elif a[0] == 1:
                    r = a[1]
                    if r == 0:
                        name = 'CM'
                    elif r in Tp:
                        name = 'C'
                    elif r in pair_of:
                        l, pp, sgn = pair_of[r]
                        name = 'B%d%d' % (y[r], y[pp])
                        coef *= sgn
                    else:
                        name = 'A%d' % y[r]
                else:
                    r = a[1]
                    name = ('DM%d' % y[0]) if r == 0 else 'D%d' % y[r]
            else:
                if a is None:
                    name = 'G%d' % y[0]
                elif a[0] == 1:
                    r = a[1]
                    if r == 0:
                        name = 'O%d' % y[0]
                    elif r in Tp:
                        name = 'C%d' % y[0]
                    elif r in pair_of:
                        l, pp, sgn = pair_of[r]
                        name = 'B%d%d%d' % (y[r], y[pp], y[0])
                        coef *= sgn
                    else:
                        name = 'A%d%d' % (y[r], y[0])
                else:
                    name = 'D%d%d' % (y[a[1]], y[0])
        vecs.setdefault(name, {})[idx] = coef
    return vecs


def form(M, va, vb):
    s = Fr(0)
    for x, cx in va.items():
        row = M[x]
        for y, cy in vb.items():
            s += cx * row[y] * cy
    return s


def perm_state(s, p, q):
    a, y = s
    sw = lambda c: q if c == p else (p if c == q else c)
    z = list(y)
    z[p], z[q] = y[q], y[p]
    return ((a[0], sw(a[1])) if a is not None else None, tuple(z))


def class_sum(st, copies, v):
    """K v for K = sum of the transpositions of the given copies; v a dict index -> coefficient"""
    idx = {s: n for n, s in enumerate(st)}
    out = {}
    for p, q in itertools.combinations(copies, 2):
        for x, c in v.items():
            y = idx[perm_state(st[x], p, q)]
            out[y] = out.get(y, 0) + c
    return {x: c for x, c in out.items() if c != 0}


def content_sum(lam):
    return sum(c - r for r, row in enumerate(lam) for c in range(row))


def dim_irrep(lam):
    lam = [x for x in lam if x > 0]
    n = sum(lam)
    if n == 0:
        return 1
    hooks = 1
    cols = [sum(1 for x in lam if x > c) for c in range(lam[0])] if lam else []
    for r, row in enumerate(lam):
        for c in range(row):
            hooks *= (row - c - 1) + (cols[c] - r - 1) + 1
    f = 1
    for x in range(2, n + 1):
        f *= x
    return f // hooks


def known_kernel(k, j, T, st, w0):
    """p (x) y for every column y of the letter-1 rows, in the rescaled frame (indicator of the column group)"""
    idx = {s: n for n, s in enumerate(st)}
    cols = sorted({y for a, y in st if a is not None and a[0] == 1})
    vecs = []
    for y in cols:
        vecs.append({idx[((1, r), y)]: Fr(1) for r in range(k)})
    return vecs
