"""Per-copy algebra for the level-2 DPS certificate.

One copy carries three qutrit factors, in the order (A, B, B') = strands 0, 1, 2.
The commutant of O -> O (x) O (x) O (O real orthogonal 3x3) on (C^3)^{(x)3} is spanned by the 15 Brauer
diagrams (perfect matchings of the 6 points: top points 0,1,2 = row indices, bottom points 3,4,5 = column
indices).  Everything here is exact (integer matrices / sympy rationals) except the *_float helpers.
"""
import itertools
from fractions import Fraction
import numpy as np
import sympy as sp

D = 3          # local dimension
NS = 3         # strands per copy: A, B, B'


def _matchings(points):
    if not points:
        yield []
        return
    a = points[0]
    for i in range(1, len(points)):
        rest = points[1:i] + points[i + 1:]
        for m in _matchings(rest):
            yield [(a, points[i])] + m


def canon(m):
    return tuple(sorted(tuple(sorted(p)) for p in m))


# fixed order: identity first, then the rest in sorted order
_ALL = sorted(set(canon(m) for m in _matchings(list(range(2 * NS)))))
IDENT = canon([(s, s + NS) for s in range(NS)])
DIAGRAMS = [IDENT] + [m for m in _ALL if m != IDENT]
NDIAG = len(DIAGRAMS)          # 15
INDEX = {m: i for i, m in enumerate(DIAGRAMS)}


def diagram_matrix(m):
    """27 x 27 integer matrix of a diagram: <i0 i1 i2| D |j0 j1 j2> = prod over pairs of deltas."""
    n = D ** NS
    M = np.zeros((n, n), dtype=np.int64)
    for idx in itertools.product(range(D), repeat=2 * NS):
        if all(idx[p] == idx[q] for p, q in m):
            r = 0
            c = 0
            for s in range(NS):
                r = r * D + idx[s]
                c = c * D + idx[NS + s]
            M[r, c] = 1
    return M


DMAT = [diagram_matrix(m) for m in DIAGRAMS]


def transpose_strands(m, strands):
    """partial transpose on the given strands: exchange top point s and bottom point s+NS."""
    sw = {}
    for s in strands:
        sw[s] = s + NS
        sw[s + NS] = s
    return canon([(sw.get(p, p), sw.get(q, q)) for p, q in m])


def pt_perm(strands):
    """permutation pi with Gamma(D_a) = D_{pi[a]}"""
    return [INDEX[transpose_strands(m, strands)] for m in DIAGRAMS]


def compose_by_matrix(a, b):
    """D_a D_b = 3^loops D_c, identified from the exact integer matrices."""
    prod = DMAT[a] @ DMAT[b]
    mx = prod.max()
    loops = 0
    while 3 ** loops < mx:
        loops += 1
    assert 3 ** loops == mx
    for c in range(NDIAG):
        if np.array_equal(prod, (3 ** loops) * DMAT[c]):
            return c, loops
    raise AssertionError('product not a diagram')


MULT = [[None] * NDIAG for _ in range(NDIAG)]   # MULT[a][b] = (c, loops)
for _a in range(NDIAG):
    for _b in range(NDIAG):
        MULT[_a][_b] = compose_by_matrix(_a, _b)


def named():
    """indices of the diagrams used below"""
    ident = INDEX[IDENT]
    F_AB = INDEX[canon([(0, 4), (1, 3), (2, 5)])]          # swap of A and B
    Om_AB = INDEX[canon([(0, 1), (3, 4), (2, 5)])]         # |Omega><Omega|_{AB} (x) I_{B'}
    S_BB = INDEX[canon([(0, 3), (1, 5), (2, 4)])]          # swap of B and B'
    F_AB2 = INDEX[canon([(0, 5), (2, 3), (1, 4)])]         # swap of A and B'
    Om_AB2 = INDEX[canon([(0, 2), (3, 5), (1, 4)])]        # |Omega><Omega|_{AB'} (x) I_B
    return dict(I=ident, F_AB=F_AB, Om_AB=Om_AB, S_BB=S_BB, F_AB2=F_AB2, Om_AB2=Om_AB2)


# ---------------------------------------------------------------- spin decomposition (exact)
def so3_generators():
    eps = np.zeros((3, 3, 3), dtype=np.int64)
    for (i, j, k), s in [((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1),
                         ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1)]:
        eps[i, j, k] = s
    L = [-eps[k] for k in range(3)]      # (L_k)_{ij} = -eps_{kij}, real antisymmetric
    return L


def total_generator(Lk, nf=NS):
    I = np.eye(D, dtype=np.int64)
    tot = np.zeros((D ** nf, D ** nf), dtype=np.int64)
    for s in range(nf):
        f = [I] * nf
        f[s] = Lk
        m = f[0]
        for x in f[1:]:
            m = np.kron(m, x)
        tot += m
    return tot


def spin_m0_bases(nf=NS):
    """exact rational bases B_j (columns) of {v : J_z v = 0, J^2 v = -j(j+1) v} (J^2 = sum_k L_k^2).
    Returns dict j -> sympy Matrix (27 x m_j)."""
    L = so3_generators()
    T = [sp.Matrix(total_generator(Lk, nf)) for Lk in L]
    Lz = T[2]
    K = T[0] ** 2 + T[1] ** 2 + T[2] ** 2
    N = sp.Matrix.hstack(*Lz.nullspace())
    # K maps ker Lz into itself: K N = N k
    k = (N.T * N).inv() * N.T * K * N
    assert K * N == N * k
    out = {}
    tot = 0
    for j in range(0, nf + 1):
        ns = (k + j * (j + 1) * sp.eye(k.shape[0])).nullspace()
        if ns:
            Bj = N * sp.Matrix.hstack(*ns)
            # clear denominators for niceness
            den = sp.ilcm(*[sp.fraction(x)[1] for x in Bj])
            Bj = Bj * den
            out[j] = Bj
            tot += Bj.shape[1]
    assert tot == N.shape[1]
    return out
