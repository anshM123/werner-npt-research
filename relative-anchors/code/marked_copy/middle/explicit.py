"""Explicit objects at small k, independent of the orbit calculus: the letter-1 block of (R_inf)_j built by j1lib.py
(Lemma B, every E^{-1} inverted sector by sector over Q), rescaled by mu (the ylib frame), the explicit twisted
vectors of mlib.NAMES as vectors over the states, and the known kernel vectors of the whole block."""
import itertools
from fractions import Fraction as Fr
from j1lib import J1, onehot, mu


def states(k, j):
    st = []
    if 0 <= j <= k:
        st += [((1, r), onehot(k, {a: 1 for a in Y})) for r in range(k) for Y in itertools.combinations(range(k), j)]
    if 0 <= j + 1 <= k:
        st += [(None, onehot(k, {a: 1 for a in Z})) for Z in itertools.combinations(range(k), j + 1)]
    return st


def block(k, j, w, bscale=1):
    """rescaled block Y~ = D (2^Z Y) D over Q at weights w (Fractions), D = diag(mu)"""
    st = states(k, j)
    J = J1(w, F=Fr, bscale=bscale)
    M = J.matrix(st)
    mus = [mu(s, w) for s in st]
    n = len(st)
    return st, [[M[a][b] * mus[a] * mus[b] for b in range(n)] for a in range(n)]


def twisted(k, j, i):
    """dict name -> vector (list over states) of the explicit twisted invariants; pairs (1,2),(3,4),.. in U"""
    st = states(k, j)
    u = k - 1
    pairs = [(1 + 2 * l, 2 + 2 * l) for l in range(i)]
    pair_of = {}
    for l, (p, q) in enumerate(pairs):
        pair_of[p] = (l, p, q, 1)
        pair_of[q] = (l, q, p, -1)
    vecs = {}
    for idx, (a, y) in enumerate(st):
        def sing(skip):
            sg = 1
            for l, (p, q) in enumerate(pairs):
                if l == skip:
                    continue
                if y[p] + y[q] != 1:
                    return 0
                sg *= 1 if y[p] == 1 else -1
            return sg
        if a is None:
            sg = sing(None)
            if sg:
                vecs.setdefault('G%d' % y[0], {})[idx] = sg
            continue
        r = a[1]
        if r == 0:
            sg = sing(None)
            if sg:
                vecs.setdefault('O%d' % y[0], {})[idx] = sg
        elif r in pair_of:
            l, rr, pp, sgn = pair_of[r]
            sg = sing(l)
            if sg:
                vecs.setdefault('B%d%d%d' % (y[rr], y[pp], y[0]), {})[idx] = sg * sgn
        else:
            sg = sing(None)
            if sg:
                vecs.setdefault('A%d%d' % (y[r], y[0]), {})[idx] = sg
    return st, vecs


def form(M, va, vb):
    s = Fr(0)
    for x, cx in va.items():
        row = M[x]
        for y, cy in vb.items():
            s += cx * row[y] * cy
    return s


def known_kernel(k, j, sig):
    """the known kernel vectors of the whole rescaled block (dicts over state indices):
    the column-group indicators; at j = k - 2 the vector sum_h s_h (vac, 1..1 - e_h); at j = k - 1 the vectors
    sum_{r,h} d_r s_h (1@r, 1..1 - e_h) for d in a basis of {d : sum_r d_r s_r = 0}"""
    st = states(k, j)
    s = [sig] + [Fr(1)] * (k - 1)
    idx = {x: n for n, x in enumerate(st)}
    vecs = []
    for Y in (itertools.combinations(range(k), j) if j >= 0 else []):
        y = onehot(k, {a: 1 for a in Y})
        vecs.append({idx[((1, r), y)]: Fr(1) for r in range(k)})
    full = tuple([1] * k)
    hole = lambda h: tuple(0 if t == h else 1 for t in range(k))
    if j == k - 2:
        vecs.append({idx[(None, hole(h))]: s[h] for h in range(k)})
    if j == k - 1:
        for r in range(1, k):
            d = [Fr(0)] * k
            d[0] = Fr(1) / s[0]
            d[r] = -Fr(1) / s[r]
            v = {}
            for rr in range(k):
                if d[rr] != 0:
                    for h in range(k):
                        v[idx[((1, rr), hole(h))]] = d[rr] * s[h]
            vecs.append(v)
    return vecs
