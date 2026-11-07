"""The irrep split of S_{U'} for the sym type I forms with a letter 1 at one marked copy (tfamlib.py), as Lemma S of
middle/NOTES.md for the letter-1 block.

The twisted vectors whose row letter is not on U' (MT: C, D_z, CM, DM_z, G; MU: O_{y0}, C_{y0}, D_{z y0}, G_{y0}) have
their U'-part in permutation modules of subsets of U' (two-row irreps of S_{U'} only).  The vectors with the row on U'
(MT: A_x, B_{xy}; MU: A_{x y0}, B_{x y y0}) form blocks (x) or (x, y0): A, B0 = B_{x 0 (y0)}, B1 = B_{x 1 (y0)}, with
n = q0 - x (MT) or n = q0 - y0 - x (MU) the Dicke count of A.  On a block with A, B0, B1 all present (a split block)
    S = i A - (n + 1) B0 + (m - n) B1,     U = A + B0,     W = B0 + B1
(the copy 0 letter y0 and the letters 2 on T' are spectators of S_{U'}, so the argument of Lemma S applies verbatim:
S spans the (u' - i - 1, i, 1) part of the block, U and W lie in the (u' - i, i) part).  Unsplit B vectors are kept.
Nothing below relies on this: the change of basis has determinant m + 1 + i on each split block (checked one-signed by
the callers), and the block-diagonality of T^T Q T is verified identically (as rational functions).  The known kernel
vectors lie in the first part: MT kappa = U_0 - W_0 + U_1 + C + CM, MU kappa_{y0} = O_{y0} + U_{0 y0} - W_{0 y0} +
U_{1 y0} + C_{y0} on split blocks (to_new checks it)."""


def blocks(case, names):
    out = []
    ys = (None,) if case == 'MT' else (0, 1)
    for x in (0, 1):
        for y0 in ys:
            if case == 'MT':
                A, B0, B1 = 'A%d' % x, 'B%d0' % x, 'B%d1' % x
            else:
                A, B0, B1 = 'A%d%d' % (x, y0), 'B%d0%d' % (x, y0), 'B%d1%d' % (x, y0)
            if all(nm in names for nm in (A, B0, B1)):
                out.append((x, y0, A, B0, B1))
    return out


def tag(x, y0):
    return '%d' % x if y0 is None else '%d%d' % (x, y0)


def basis(case, names, q0, m, i, one):
    """new basis vectors as dicts old name -> coefficient; returns (vecs, two-part names, three-part names, split)"""
    vecs, used = {}, set()
    split = blocks(case, names)
    three, newtwo = [], []
    for x, y0, A, B0, B1 in split:
        n = q0 - (x if y0 is None else x + y0)
        tg = tag(x, y0)
        vecs['U' + tg] = {A: one, B0: one}
        vecs['W' + tg] = {B0: one, B1: one}
        vecs['S' + tg] = {A: i * one, B0: -(n + 1) * one, B1: (m - n) * one}
        newtwo += ['U' + tg, 'W' + tg]
        three.append('S' + tg)
        used |= {A, B0, B1}
    two = []
    for nm in names:
        if nm not in used:
            vecs[nm] = {nm: one}
            two.append(nm)
    two += newtwo
    return vecs, two, three, split


def congruence(names, Q, vecs, order, zero):
    """T^T Q T on the new basis, in the given order of new names: a dict (p, q) -> entry"""
    idx = {nm: a for a, nm in enumerate(names)}
    QT = {}
    for nb in order:
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
    for x, na in enumerate(order):
        for nb in order[x:]:
            s = zero
            for oa, cf in vecs[na].items():
                v = QT[nb][idx[oa]]
                if v != 0:
                    s = s + cf * v
            M[na, nb] = s
            M[nb, na] = s
    return M


def to_new(v, split):
    """a vector of the first part, given in the old basis (old name -> coefficient), in the new basis; on a split block
    (a, b, c) on (A, B0, B1) becomes a U + c W, which requires b = a + c (None otherwise)"""
    out, done = {}, set()
    for x, y0, A, B0, B1 in split:
        a, b, c = v.get(A, 0), v.get(B0, 0), v.get(B1, 0)
        if (a + c) - b != 0:
            return None
        tg = tag(x, y0)
        if a != 0:
            out['U' + tg] = a
        if c != 0:
            out['W' + tg] = c
        done |= {A, B0, B1}
    for nm, cf in v.items():
        if nm not in done:
            out[nm] = cf
    return out
