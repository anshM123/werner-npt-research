"""Torus blocks of (R_inf)_j for any k, in the rescaled basis of ylib.py, with a gauge fixing and exact tools.

Block label: (n1, n2_0, ..., n2_{k-1}) = wstate(a, y) (total number of letters 1, and the number of letters 2 at each
copy, row and column together).  Kinds: 'L1' (no letter 2), 'I' (every n2_t <= 1, some = 1), 'II' (some n2_t = 2).

Gauge: the kernel of the block is computed exactly at a rational point; coordinates are deleted greedily in a fixed
preference order ((1@g, .) states first with g = 1, then the other row letters 1, then vacuum states, then the rest)
so that the kernel restricted to the deleted set is invertible.  For the known kernels (README, Section 3; ../middle_k4/NOTES.md,
Section 2) this reproduces the gauge of k4lib.py there.  The gauge-fixed matrix is the principal submatrix
on the kept coordinates; it is positive definite iff the block is >= 0 with kernel exactly the computed kernel.
"""
import itertools, random
import flint
from ylib import Builder, Ring, rinf_states, wstate, point_ring
from flintfield import Field

Q = flint.fmpq


def blocks(k, j):
    """dict label -> list of states (a, y), and kind per label."""
    out = {}
    for st in rinf_states(k, j):
        out.setdefault(wstate(*st), []).append(st)
    kinds = {}
    for lab in out:
        n2 = lab[1:]
        if all(x == 0 for x in n2):
            kinds[lab] = 'L1'
        elif any(x == 2 for x in n2):
            kinds[lab] = 'II'
        else:
            kinds[lab] = 'I'
    return out, kinds


def canon_T(lab):
    return tuple(t for t, x in enumerate(lab[1:]) if x == 1)


def representatives(k, j, kinds_wanted=('L1', 'I')):
    """one block per orbit under permutations of the copies (orbit of the type I label is fixed by |T|; for the
    letter-1 block there is one).  Type II blocks: representative per pattern of n2 up to permutation."""
    B, kinds = blocks(k, j)
    seen = set()
    reps = []
    for lab in sorted(B):
        if kinds[lab] not in kinds_wanted:
            continue
        key = (lab[0], tuple(sorted(lab[1:])))
        if key in seen:
            continue
        # prefer the label whose letter-2 copies are the first ones (T = {0, 1, ...})
        n2s = sorted(lab[1:], reverse=True)
        lab2 = (lab[0],) + tuple(n2s)
        if lab2 in B:
            lab = lab2
        seen.add(key)
        reps.append((lab, kinds[lab]))
    return reps, B


def label_state(st):
    a, y = st
    s2 = lambda z: ''.join(str(c) for c in z)
    return f'({s2(a)},{s2(y)})'


def sym_builder(k, dehom=False):
    if dehom:
        K = Field([f's{t}' for t in range(1, k)])
        R = Ring([K(1)] + list(K.gens), K(0), K(1))
        s = [K(1)] + list(K.gens)
    else:
        K = Field([f's{t}' for t in range(k)])
        R = Ring(list(K.gens), K(0), K(1))
        s = list(K.gens)
    return Builder(R), K, s


def pt_builder(svals):
    return Builder(point_ring([Q(x) for x in svals]))


def nullspace_Q(Mx):
    """exact nullspace basis (list of vectors) of a square fmpq matrix (list of lists)."""
    n = len(Mx)
    A = flint.fmpq_mat(n, n, [Mx[i][j] for i in range(n) for j in range(n)])
    R, rank = A.rref()
    piv = []
    r = 0
    for c in range(n):
        if r < rank and R[r, c] != 0:
            piv.append(c)
            r += 1
    free = [c for c in range(n) if c not in piv]
    basis = []
    for f in free:
        v = [Q(0)] * n
        v[f] = Q(1)
        for i, pc in enumerate(piv):
            v[pc] = -R[i, f]
        basis.append(v)
    return basis


def rank_Q(rows):
    if not rows:
        return 0
    m, n = len(rows), len(rows[0])
    A = flint.fmpq_mat(m, n, [x for r in rows for x in r])
    return A.rank()


def gauge(states, Mpt, g=1):
    """greedy gauge: returns (kernel basis at the point, deleted indices, kept indices)."""
    ker = nullspace_Q(Mpt)
    n = len(states)
    k = len(states[0][0])
    vac = (0,) * k

    def pref(i):
        a, y = states[i]
        if a != vac and a[g] == 1:
            return (0, i)
        if a != vac and 1 in a:
            return (1, a.index(1), i)
        if a == vac:
            return (2, i)
        return (3, i)
    order = sorted(range(n), key=pref)
    dele = []
    for i in order:
        if len(dele) == len(ker):
            break
        cand = dele + [i]
        if rank_Q([[v[c] for c in cand] for v in ker]) == len(cand):
            dele = cand
    assert len(dele) == len(ker)
    keep = [i for i in range(n) if i not in dele]
    return ker, dele, keep


def rand_point(k, rng, lo=1, hi=40):
    return [Q(rng.randint(lo, hi), rng.randint(lo, hi)) for _ in range(k)]
