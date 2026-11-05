"""Torus blocks of (R_inf)_1 and (R_inf)_2 at k = 4, built with ylib.Builder (copied from relative-anchors/code/).

ylib works in the rescaled basis x = Psi / mu, mu_X = prod_t w_t^{N_t(X)} (N_t = number of letters 1 at copy t in row
and column of the state X), where every entry is a rational function of the squared weights s_t = w_t^2.

A torus block is fixed by j and by the set T of copies carrying one letter 2 (T = () is the letter-1 block).  The
total number of letters 1 in the block is n1 = j + 1 - |T|, the same for all its states; under s -> lam s every entry
of the rescaled block is multiplied by lam^{n1}, so entries are homogeneous of degree n1 in s and we may put s0 = 1.

Known kernel vectors, in rescaled coordinates:
  p (x) y  (y a column string of the block)  ->  the all-ones vector on the states (1@r, y), r = 0..3;
  vac (x) J^- |1111>  (j = 2, letter-1 block)  ->  x_(vac, hole h) = s_h  (zero elsewhere).
Gauge: delete (1@g, y) for one fixed g and every y, and for j = 2, T = () also the state (vac, hole g0) with s_g0 != 0.
"""
import itertools
import flint
from ylib import Builder, Ring, rinf_states, wstate, point_ring
from flintfield import Field

K4 = 4


def block_states(j, T, k=K4):
    T = set(T)
    out = []
    for (a, y) in rinf_states(k, j):
        n2 = wstate(a, y)[1:]
        if all(n2[t] == (1 if t in T else 0) for t in range(k)):
            out.append((a, y))
    return out


def label(st):
    a, y = st
    def s2(z):
        return ''.join(str(c) for c in z)
    return f'({s2(a)},{s2(y)})'


def kernel_vectors(states, j, T, s):
    """the known kernel vectors (rescaled coordinates) as lists over states; s = squared weights (ring elements)."""
    k = len(states[0][0])
    vac = (0,) * k
    vecs = []
    cols = sorted({y for (a, y) in states if a != vac and sum(a) == 1 and 1 in a})
    for y in cols:
        v = [1 if (a != vac and 1 in a and y2 == y) else 0 for (a, y2) in states]
        if sum(v) == k:
            vecs.append(v)
    if j == k - 2 and len(T) == 0:
        v = []
        for (a, y) in states:
            if a == vac and y.count(0) == 1:
                v.append(s[y.index(0)])
            else:
                v.append(0)
        vecs.append(v)
    return vecs


def gauge_keep(states, j, T, g=1, g0=0):
    """indices kept after deleting the gauge coordinates (1@g, y) for every y, and (vac, hole g0) if j = k-2, T = ()."""
    k = len(states[0][0])
    vac = (0,) * k
    keep = []
    for i, (a, y) in enumerate(states):
        if a != vac and a[g] == 1:
            continue
        if j == k - 2 and len(T) == 0 and a == vac and y.count(0) == 1 and y.index(0) == g0:
            continue
        keep.append(i)
    return keep


def sym_builder(dehom=True):
    """Builder over Q(s1, s2, s3) with s0 = 1 (dehom) or over Q(s0, ..., s3)."""
    if dehom:
        K = Field(['s1', 's2', 's3'])
        R = Ring([K(1)] + list(K.gens), K(0), K(1))
    else:
        K = Field(['s0', 's1', 's2', 's3'])
        R = Ring(list(K.gens), K(0), K(1))
    return Builder(R), K


def pt_builder(svals):
    return Builder(point_ring(svals))


BLOCKS = [(1, ()), (1, (0,)), (1, (0, 1)), (2, ()), (2, (0,)), (2, (0, 1)), (2, (0, 1, 2))]
