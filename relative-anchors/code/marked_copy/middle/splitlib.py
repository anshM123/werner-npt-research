"""Splitting the twisted space V_i into the parts of the two S_u-irreps (u - i, i) and (u - i - 1, i, 1).

For x, y0 in {0, 1} put n = q0 - y0 - x.  On span(A_{x,y0}, B_{x,0,y0}, B_{x,1,y0}) (row letter in U) the equivariant
maps (1@r, Y) -> Y - r and (1@r, Y) -> Y + r (into permutation modules of subsets, which carry only two-row irreps)
send A, B0, B1 to multiples of singlets x Dicke: phi1 = ((m - n), 0, -i), phi2 = ((n + 1), i, 0).  When A_{x,y0} is
present (0 <= n <= m - 1, and then B0, B1 are present for i >= 1 and both images are nonzero) the common kernel is
spanned by
    S_{x,y0} = i A_{x,y0} - (n + 1) B_{x,0,y0} + (m - n) B_{x,1,y0},
and its orthogonal complement (for the norms |A|^2 : |B0|^2 : |B1|^2 = (m - n)/i : (m - n)/(n + 1) : 1) is spanned by
    U_{x,y0} = A_{x,y0} + B_{x,0,y0},     W_{x,y0} = B_{x,0,y0} + B_{x,1,y0}
(such a block (x, y0) is called split).  When A_{x,y0} is absent, the B vectors present are kept as they are, in the
first part.  The new basis is O.., G.., U.., W.. and the unsplit B.. (the (u - i, i) part) and S.. (the
(u - i - 1, i, 1) part).  Nothing below relies on representation theory: on a split block the change of basis has
determinant m + 1 + i (> 0), elsewhere it is the identity, and the block-diagonality of T^T Q T is verified identically
(as rational functions) by the caller."""

TWO_ORDER = ['O0', 'O1', 'G0', 'G1', 'U00', 'U01', 'U10', 'U11', 'W00', 'W01', 'W10', 'W11',
             'B000', 'B001', 'B010', 'B011', 'B100', 'B101', 'B110', 'B111', 'A00', 'A01', 'A10', 'A11']
THREE = ['S00', 'S01', 'S10', 'S11']


def blocks(names):
    """the split blocks (x, y0): A_{x,y0}, B_{x,0,y0}, B_{x,1,y0} all present"""
    out = []
    for x in (0, 1):
        for y0 in (0, 1):
            if all(nm in names for nm in ('A%d%d' % (x, y0), 'B%d0%d' % (x, y0), 'B%d1%d' % (x, y0))):
                out.append((x, y0))
    return out


def basis(names, q0, m, i, one):
    """the new basis vectors as dicts old name -> coefficient (field elements built from q0, m, i, one)"""
    vecs = {}
    split = blocks(names)
    used = set()
    for x, y0 in split:
        n = q0 - (y0 + x)
        A, B0, B1 = 'A%d%d' % (x, y0), 'B%d0%d' % (x, y0), 'B%d1%d' % (x, y0)
        vecs['U%d%d' % (x, y0)] = {A: one, B0: one}
        vecs['W%d%d' % (x, y0)] = {B0: one, B1: one}
        vecs['S%d%d' % (x, y0)] = {A: i * one, B0: -(n + 1) * one, B1: (m - n) * one}
        used |= {A, B0, B1}
    for nm in names:
        if nm not in used:
            vecs[nm] = {nm: one}
    return vecs, split


def new_names(vecs):
    two = [nm for nm in TWO_ORDER if nm in vecs]
    three = [nm for nm in THREE if nm in vecs]
    assert len(two) + len(three) == len(vecs)
    return two, three


def congruence(names, Q, vecs, zero):
    """T^T Q T for the basis vecs (dict new name -> dict old name -> coefficient): returns (new names, matrix)"""
    two, three = new_names(vecs)
    new = two + three
    idx = {nm: a for a, nm in enumerate(names)}
    QT = {}
    for nb in new:
        col = []
        for a in range(len(names)):
            s = zero
            for ob, cf in vecs[nb].items():
                q = Q[a][idx[ob]]
                if q != 0:
                    s = s + q * cf
            col.append(s)
        QT[nb] = col
    M = {}
    for x, na in enumerate(new):
        for nb in new[x:]:
            s = zero
            for oa, cf in vecs[na].items():
                v = QT[nb][idx[oa]]
                if v != 0:
                    s = s + cf * v
            M[na, nb] = s
            M[nb, na] = s
    return two, three, M


def to_new(v, split):
    """express a vector of the first part, given in the old basis (dict old name -> coefficient), in the new basis:
    on a split block (a, b, c) on (A, B0, B1) becomes a U + c W, which requires b = a + c (checked)"""
    out = {}
    done = set()
    for x, y0 in split:
        A, B0, B1 = 'A%d%d' % (x, y0), 'B%d0%d' % (x, y0), 'B%d1%d' % (x, y0)
        a, b, c = v.get(A, 0), v.get(B0, 0), v.get(B1, 0)
        if (a + c) - b != 0:
            return None
        if a != 0:
            out['U%d%d' % (x, y0)] = a
        if c != 0:
            out['W%d%d' % (x, y0)] = c
        done |= {A, B0, B1}
    for nm, cf in v.items():
        if nm not in done:
            out[nm] = cf
    return out
