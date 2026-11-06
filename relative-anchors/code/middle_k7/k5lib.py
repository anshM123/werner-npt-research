"""Blocks of (R_inf)_j at k copies with two weight values: states, known kernel, gauge, configurations, builders.

Families (k = 5 unless stated):
  L3   letter-1 block of (R_inf)_3 (T = {});      known kernel p (x) |1@Y> (|Y| = 3) and vac (x) J^-|11111> (K2)
  T31  type I block of (R_inf)_3 with |T| = 1;   known kernel p (x) y for every column y of the letter-1 rows
  T32  type I block of (R_inf)_3 with |T| = 2;   same
  T33  type I block of (R_inf)_3 with |T| = 3;   same
  T42  type I block of (R_inf)_4 with |T| = 2;   same (p (x) f_h, h in U)
General: family 'j,m' means the type I block of (R_inf)_j with |T| = m (m = 0: the letter-1 block).

States are ylib pairs (a, y): a the row string (all zeros for the vacuum row), y the column string, letters 0, 1, 2.
The rescaled frame is that of ylib.py: Yr = D Y D with D = diag(mu), mu = prod_t w_t^{N_t} (N_t = number of letters 1
at copy t in row and column).  ylib.py builds Yr without the factor 2^Z of Lemma B (Z = k - j).

Configurations.  The block labelled by T is mapped by a permutation of the copies to the block labelled by
T0 = {0, ..., m-1}; U0 = {m, ..., k-1}.  A two-value configuration is (a, b): the copies 0..a-1 of T0 and m..m+b-1 of
U0 have squared weight 1, the others squared weight q.  (For m = 0, b is the number of copies of weight 1.)
"""
import itertools
from fractions import Fraction as Fr
import flint
import ylib
from j1lib import J1

FAMILIES = {'L3': (3, 0), 'T31': (3, 1), 'T32': (3, 2), 'T33': (3, 3), 'T42': (4, 2)}


def fam(name):
    """(j, m) of a family: the names above, or 'jJmM' (e.g. j3m0: the letter-1 block of (R_inf)_3)"""
    if name in FAMILIES:
        return FAMILIES[name]
    if name.startswith('j') and 'm' in name:
        j, m = name[1:].split('m')
        return int(j), int(m)
    j, m = name.split(',')
    return int(j), int(m)


def block_states(k, j, T):
    """states of the torus block of (R_inf)_j labelled by T (each copy of T carries the letter 2 exactly once, in row
    or column; no other letter 2).  Order: letter-2 rows (by copy, then column), letter-1 rows (by copy, then
    column), vacuum states (by column)."""
    T = tuple(sorted(T))
    target = (j + 1 - len(T),) + tuple(1 if t in T else 0 for t in range(k))
    sts = [st for st in ylib.rinf_states(k, j) if ylib.wstate(*st) == target]

    def key(st):
        a, y = st
        if any(a):
            r = next(i for i in range(k) if a[i])
            return (0 if a[r] == 2 else 1, r, y)
        return (2, 0, y)
    return sorted(sts, key=key)


def row_of(st):
    a, y = st
    for i, x in enumerate(a):
        if x:
            return (x, i)
    return None


def kernel_groups(k, states):
    """columns y that carry the letter-1 rows: for each, the indices of the states (1@r, y), r = 0..k-1"""
    groups = {}
    for i, st in enumerate(states):
        rl = row_of(st)
        if rl is not None and rl[0] == 1:
            groups.setdefault(st[1], []).append(i)
    for y, idx in groups.items():
        assert len(idx) == k, (y, idx)
    return groups


def known_kernel(k, j, T, states, s, zero, one, r0=None):
    """known kernel vectors in the rescaled frame and the deleted coordinate of each.
    p (x) y  ->  indicator of the column group of y (deleted: the state (1@r0, y));
    for the letter-1 block of (R_inf)_{k-2}: vac (x) J^-|1...1>  ->  s_h on (vac, 1...1 - e_h) (deleted: h = 0)."""
    r0 = k - 1 if r0 is None else r0
    n = len(states)
    vecs, dels = [], []
    groups = kernel_groups(k, states)
    for y in sorted(groups):
        v = [zero] * n
        for i in groups[y]:
            v[i] = one
        vecs.append(v)
        d = [i for i in groups[y] if row_of(states[i])[1] == r0]
        assert len(d) == 1
        dels.append(d[0])
    if len(T) == 0 and j == k - 2:
        v = [zero] * n
        idx = {}
        for i, (a, y) in enumerate(states):
            if not any(a) and sum(1 for x in y if x == 1) == k - 1:
                h = next(t for t in range(k) if y[t] == 0)
                idx[h] = i
        assert len(idx) == k
        for h, i in idx.items():
            v[i] = s[h]
        vecs.append(v)
        dels.append(idx[0])
    # each known vector is the only one nonzero on its deleted coordinate
    for a_, d in enumerate(dels):
        for b_, v in enumerate(vecs):
            assert (v[d] != 0) == (a_ == b_)
    return vecs, dels


def config_s(k, m, a, b, qv, one):
    """squared weights of configuration (a, b) at ratio qv (element of the field)"""
    s = []
    for t in range(k):
        if t < m:
            s.append(one if t < a else qv)
        else:
            s.append(one if t - m < b else qv)
    return s


def configs(k, m):
    """all configurations (a, b), 0 <= a <= m, 0 <= b <= k - m"""
    return [(a, b) for a in range(m + 1) for b in range(k - m + 1)]


def canonical(k, m, a, b):
    """the representative of (a, b) under the level exchange (a, b) ~ (m - a, k - m - b) (at ratio 1/q)"""
    return min((a, b), (m - a, k - m - b))


def rep_configs(k, m):
    return sorted(set(canonical(k, m, a, b) for a, b in configs(k, m)))


# ---------------------------------------------------------------- builders
def build_ylib(k, j, T, s, zero, one, bscale=1, states=None):
    """the rescaled block (ylib units) over the ring of s; bscale multiplies the B-part"""
    states = block_states(k, j, T) if states is None else states
    R = ylib.Ring(s, zero, one)
    bld = ylib.Builder(R)
    n = len(states)
    M = [[None] * n for _ in range(n)]
    for i in range(n):
        a, yp = states[i]
        for l in range(i, n):
            b, y = states[l]
            v = bld.A(a, yp, b, y)
            bb = bld.B(a, yp, b, y)
            if bscale == 1:
                v = v - bb
            else:
                v = v - bb * bscale
            M[i][l] = v
    for i in range(n):
        for l in range(i):
            M[i][l] = M[l][i]
    return states, M


def to_j1(st, k):
    a, y = st
    rl = row_of(st)
    return (rl, y)


def mu(st, w):
    a, y = st
    m = Fr(1)
    for t in range(len(y)):
        m *= w[t] ** ((1 if a[t] == 1 else 0) + (1 if y[t] == 1 else 0))
    return m


def build_j1(k, j, T, w, bscale=1, states=None):
    """the rescaled block in ylib units from Lemma B (j1lib.py) at rational weights w:
    Yr[X1, X2] = mu(X1) mu(X2) (2^Z Y)[X1, X2] / 2^Z"""
    states = block_states(k, j, T) if states is None else states
    J = J1([Fr(x) for x in w], bscale=Fr(bscale))
    Z = k - j
    n = len(states)
    M = [[None] * n for _ in range(n)]
    for i in range(n):
        for l in range(i, n):
            v = J.entry(to_j1(states[i], k), to_j1(states[l], k)) * mu(states[i], w) * mu(states[l], w) / 2 ** Z
            M[i][l] = M[l][i] = Fr(v)
    return states, M


def rsqrt(x):
    """the rational square root of a rational square"""
    from math import isqrt
    x = Fr(x)
    n, d = isqrt(x.numerator), isqrt(x.denominator)
    assert n * n == x.numerator and d * d == x.denominator, x
    return Fr(n, d)


# ---------------------------------------------------------------- exact linear algebra over Q
def to_fmpq_mat(M):
    n = len(M)
    return flint.fmpq_mat(n, len(M[0]), [flint.fmpq(Fr(x).numerator, Fr(x).denominator) for row in M for x in row])


def rank_q(M):
    return to_fmpq_mat(M).rank()


def ldl_pd(M):
    """exact: True if the symmetric rational matrix M is positive definite (all pivots of symmetric elimination > 0);
    returns (ok, index of the first nonpositive pivot or None)"""
    A = [[Fr(x) for x in row] for row in M]
    n = len(A)
    for c in range(n):
        p = A[c][c]
        if p <= 0:
            return False, c
        for i in range(c + 1, n):
            if A[i][c] != 0:
                f = A[i][c] / p
                for l in range(c + 1, n):
                    if A[c][l] != 0:
                        A[i][l] -= f * A[c][l]
    return True, None


def matvec(M, v):
    return [sum((M[i][l] * v[l] for l in range(len(v)) if v[l] != 0), 0) for i in range(len(M))]


def principal(M, keep):
    return [[M[i][l] for l in keep] for i in keep]
