"""Reduced, gauge-fixed forms of the letter-1 block of (R_inf)_2 at fixed level sizes (u0 copies of squared weight 1,
u1 copies of squared weight q), one per irrep of S_u0 x S_u1, over any field of numbers (Fractions at a point, or the
rational functions of ufield.RF).

For an irrep lam with Young symmetrizer e (young.py): KP = basis of e.K (kernel span), WP = basis of e.W (a complement
of e.K in e.M_P), QB = basis of e.M_Q.  With Y_QQ = 2 What3^{-1} (Lemma G, b = 1) the Schur complement of the vacuum
states restricted to the part is
    S = WP^T Y_PP WP - 1/2 R_QW^T G_Q^{-1} (QB^T What3 QB) G_Q^{-1} R_QW,     R_QW = QB^T Y_QP WP,  G_Q = QB^T QB,
and for b != 1 (negative controls) S = WP^T Y_PP WP - R_QW^T R_QQ^{-1} R_QW with
    R_QQ = 2(1 - b) QB^T D2 QB + 2b G_Q (QB^T What3 QB)^{-1} G_Q   (D2 = diag(pi_T)).
The letter-1 block is >= 0 with kernel exactly span{p (x) |1@P>} iff, for every lam, the kernel identities hold and
S is positive definite."""
import itertools
from fractions import Fraction as Fr
from l2block import Block, ginv, msets
from young import part_data, irreps


def matmul(A, B, zero):
    n, m, p = len(A), len(B), len(B[0]) if B else 0
    out = [[zero] * p for _ in range(n)]
    for i in range(n):
        Ai = A[i]
        for l in range(m):
            a = Ai[l]
            if a == 0:
                continue
            Bl = B[l]
            row = out[i]
            for j in range(p):
                if not (Bl[j] == 0):
                    row[j] = row[j] + a * Bl[j]
    return out


def transpose(A):
    return [list(r) for r in zip(*A)] if A else []


def apply_sparse(M, vecs, zero):
    """M (dense n x n over the field) times each integer vector of vecs (dense lists); returns n x len(vecs)"""
    sp = [[(j, c) for j, c in enumerate(v) if c != 0] for v in vecs]
    out = []
    for row in M:
        out.append([sum((row[j] * c for j, c in s), zero) for s in sp])
    return out


def tdot(vecs, cols, zero):
    """vecs (integer dense lists) transposed times cols (n x m matrix over the field)"""
    out = []
    for v in vecs:
        sp = [(j, c) for j, c in enumerate(v) if c != 0]
        out.append([sum((cols[j][l] * c for j, c in sp), zero) for l in range(len(cols[0]) if cols else 0)])
    return out


def levels_s(u0, u1, q, one):
    return [one] * u0 + [q * one] * u1


def build_block(u0, u1, q, one, b=1):
    return Block(levels_s(u0, u1, q, one), one, bscale=b, need_qq=False)


def reduce_all(blk, u0, u1, b=1, log=None, parts=None):
    """returns (list of per-irrep dicts, kernel_ok) for the block object blk (l2block.Block, need_qq not needed)"""
    k = blk.k
    one, zero = blk.one, blk.zero
    PP = blk.PP()
    PQ = blk.PQ()
    W3 = blk.What3()
    QP = transpose(PQ)
    nP, nQ = len(PP), len(QP)
    # kernel identities on the whole block: Y_PP 1_X = 0 and Y_QP 1_X = 0 for every pair X
    kern_ok = True
    if b == 1:
        for X in msets(k, 2):
            cols = [i for i, (r, XX) in enumerate(blk.P) if XX == X]
            for i in range(nP):
                if not (sum((PP[i][j] for j in cols), zero) == 0):
                    kern_ok = False
            for i in range(nQ):
                if not (sum((QP[i][j] for j in cols), zero) == 0):
                    kern_ok = False
    D2 = [blk.pi(T) for T in blk.Qs]
    out = []
    for (l0, l1) in irreps(u0, u1):
        if parts is not None and (l0, l1) not in parts:
            continue
        pd = part_data(k, u0, l0, l1)
        if pd['mP'] == 0 and len(pd['Q']) == 0:
            continue
        rec = dict(lam=(l0, l1), dim=pd['dim'], mP=pd['mP'], d=pd['d'], nW=len(pd['W']), nQ=len(pd['Q']))
        WP, KP, QB = pd['W'], pd['K'], pd['Q']
        # reduced kernel checks (implied by the full identities, re-checked in reduced form)
        if b == 1 and KP:
            YK = apply_sparse(PP, KP, zero)
            rec['ker_red_ok'] = all(x == 0 for row in tdot(KP + WP, YK, zero) for x in row)
        if not WP:
            rec['S'] = []
            out.append(rec)
            continue
        YW = apply_sparse(PP, WP, zero)
        RWW = tdot(WP, YW, zero)
        if QB:
            QPW = apply_sparse(QP, WP, zero)              # nQ x nW
            RQW = tdot(QB, QPW, zero)                      # mQ x nW
            GQ = [[sum((Fr(a) * c for a, c in zip(u, v)), Fr(0)) for v in QB] for u in QB]
            W3Q = apply_sparse(W3, QB, zero)
            RWh = tdot(QB, W3Q, zero)
            GQi = [[one * x for x in row] for row in ginv([[Fr(x) for x in r] for r in GQ], Fr(1))]
            if b == 1:
                mid = matmul(matmul(GQi, RWh, zero), GQi, zero)
                corr = matmul(matmul(transpose(RQW), mid, zero), RQW, zero)
                S = [[RWW[i][j] - corr[i][j] / 2 for j in range(len(WP))] for i in range(len(WP))]
            else:
                GQf = [[one * x for x in row] for row in GQ]
                D2Q = [[sum((D2[t] * (QB[i][t] * QB[j][t]) for t in range(nQ) if QB[i][t] != 0 and QB[j][t] != 0),
                            zero) for j in range(len(QB))] for i in range(len(QB))]
                RWhi = ginv(RWh, one)
                RQQ = [[2 * (1 - b) * D2Q[i][j] for j in range(len(QB))] for i in range(len(QB))]
                t2 = matmul(matmul(GQf, RWhi, zero), GQf, zero)
                RQQ = [[RQQ[i][j] + 2 * b * t2[i][j] for j in range(len(QB))] for i in range(len(QB))]
                corr = matmul(matmul(transpose(RQW), ginv(RQQ, one), zero), RQW, zero)
                S = [[RWW[i][j] - corr[i][j] for j in range(len(WP))] for i in range(len(WP))]
        else:
            S = RWW
        rec['S'] = S
        out.append(rec)
        if log:
            log(f'  part {rec["lam"]}: dim {rec["dim"]}, mult in M_P {rec["mP"]}, in K {rec["d"]}, gauge-fixed size '
                f'{rec["nW"]}, Q-vectors {rec["nQ"]}')
    return out, kern_ok


def leading_pivots(S, one):
    """symmetric elimination without pivoting; returns the pivots (None where a pivot vanishes)"""
    A = [list(r) for r in S]
    n = len(A)
    piv = []
    for c in range(n):
        p = A[c][c]
        if p == 0:
            piv.append(None)
            return piv
        piv.append(p)
        for i in range(c + 1, n):
            if not (A[i][c] == 0):
                f = A[i][c] / p
                for j in range(c + 1, n):
                    if not (A[c][j] == 0):
                        A[i][j] = A[i][j] - f * A[c][j]
    return piv
