"""Exact validation of the reduced forms of l1sym.py against explicit projections of the letter-1 block of (R_inf)_1
built from Lemma B by j1lib.py (rescaled frame), at rational weights and integer level sizes.

For each point (w0, w1, u0, u1) and bscale b in {1, 1001/1000} and each part that exists at these sizes:
 (a) the explicit test vectors are built on the concrete copies (named copies of the seed = the first copies of their
     level) and R_PP, R_PQ, R_QQ = V^T Y~ V are computed from the j1lib block; they must equal the pattern sums of
     l1sym.py evaluated in exact rational arithmetic;
 (b) the Schur complement of the vacuum states of the j1lib block, projected on the P-vectors, must equal S;
 (c) the Gram matrix of the kept vectors must be nonsingular (the vectors are a basis of the multiplicity space);
 (d) the multiplicities times the irrep dimensions must add up to k^2 (P-states) and k(k-1)/2 (Q-states);
 (e) the kernel vectors p (x) |1@c> of the block are annihilated by S (b = 1).
Usage: python check_reduce.py"""
import sys, time, itertools
from fractions import Fraction as Fr
from j1lib import J1, block_states, mu
from l1sym import Model, part_spec, reduce_part, keep_for, mm, mT, minv, det

t0 = time.time()


def concrete_labels(part, u0, u1, swap=False):
    """labels (level, name, copy) of the concrete copies; level 0 first, or (swap) level 1 first."""
    named = part_spec(part)[0]
    k = u0 + u1
    levs = ([1] * u1 + [0] * u0) if swap else ([0] * u0 + [1] * u1)
    names = [None] * k
    nxt = {0: u1, 1: 0} if swap else {0: 0, 1: u0}
    for (lev, nm, _) in named:
        names[nxt[lev]] = nm
        nxt[lev] += 1
    return [(levs[t], names[t], t) for t in range(k)]


def explicit(part, w0, w1, u0, u1, b, swap=False):
    k = u0 + u1
    w = ([w1] * u1 + [w0] * u0) if swap else ([w0] * u0 + [w1] * u1)
    J = J1(w, bscale=b)
    st = block_states(k, ())
    M = J.matrix(st)
    m = [mu(z, J.w) for z in st]
    Mt = [[m[i] * M[i][j] * m[j] for j in range(len(st))] for i in range(len(st))]
    labs = concrete_labels(part, u0, u1, swap)
    named, P, Qv, ker, names = part_spec(part)
    kp, kq, gauge = keep_for(part, u0, u1)
    P = [P[i] for i in kp]
    Qv = [Qv[i] for i in kq] if (Qv and kq is not None) else []
    Pi = [i for i, (a, y) in enumerate(st) if a is not None]
    Qi = [i for i, (a, y) in enumerate(st) if a is None]

    def vecP(f):
        v = [Fr(0)] * len(st)
        for i in Pi:
            a, y = st[i]
            v[i] = Fr(f(labs[a[1]], labs[y.index(1)]))
        return v

    def vecQ(g):
        v = [Fr(0)] * len(st)
        for i in Qi:
            pr = [t for t in range(k) if st[i][1][t] == 1]
            v[i] = Fr(g(labs[pr[0]], labs[pr[1]]))
        return v
    VP = [vecP(f) for f in P]
    VQ = [vecQ(g) for g in Qv]

    def form(X, Y, A):
        return [[sum(x[i] * A[i][j] * y[j] for i in range(len(st)) if x[i] for j in range(len(st)) if y[j])
                 for y in Y] for x in X]
    out = dict(RPP=form(VP, VP, Mt), GP=[[sum(a * c for a, c in zip(x, y)) for y in VP] for x in VP])
    if VQ:
        out['RPQ'] = form(VP, VQ, Mt)
        out['RQQ'] = form(VQ, VQ, Mt)
    # Schur complement of the vacuum states, projected
    nQ = len(Qi)
    A = [[Mt[i][j] for j in Qi] for i in Qi]
    Ai = minv(A, Fr(1))
    BPQ = [[Mt[i][j] for j in Qi] for i in Pi]
    C = mm(mm(BPQ, Ai, Fr(0)), mT(BPQ), Fr(0))
    SP = [[Mt[Pi[a]][Pi[c]] - C[a][c] for c in range(len(Pi))] for a in range(len(Pi))]
    VPr = [[v[i] for i in Pi] for v in VP]
    out['S'] = [[sum(x[a] * SP[a][c] * y[c] for a in range(len(Pi)) if x[a] for c in range(len(Pi)) if y[c])
                 for y in VPr] for x in VPr]
    return out


def dims(u):
    """dimensions of the irreps (u), (u-1,1), (u-2,2), (u-2,1,1) of S_u (0 if absent)."""
    return dict(inv=1, std=u - 1 if u >= 2 else 0, two=u * (u - 3) // 2 if u >= 4 else 0,
                anti=(u - 1) * (u - 2) // 2 if u >= 3 else 0)


def count_check(u0, u1):
    """sum over parts of dim(irrep) * multiplicity equals k^2 (P) and k(k-1)/2 (Q)."""
    k = u0 + u1
    totP = totQ = 0
    d0, d1 = dims(u0), dims(u1)
    for part in ('inv', 'std', 'two', 'anti', 'mix'):
        for lev in (0, 1):
            if part == 'inv' and lev == 1:
                continue
            if part == 'mix' and lev == 1:
                continue
            ua, ub = (u0, u1) if lev == 0 else (u1, u0)
            da = dims(ua)
            if part == 'mix':
                dim = d0['std'] * d1['std']
            else:
                dim = da[part]
            if dim == 0:
                continue
            kp, kq, g = keep_for(part, ua, ub)
            nq = len(kq) if kq else 0
            totP += dim * len(kp)
            totQ += dim * nq
    return totP == k * k and totQ == k * (k - 1) // 2, (totP, totQ)


def main():
    ok = True
    pts = [(Fr(1), Fr(1, 3), 2, 2), (Fr(2), Fr(1), 3, 1), (Fr(1), Fr(5, 2), 1, 3), (Fr(1, 2), Fr(3), 3, 2),
           (Fr(3), Fr(1, 4), 2, 3), (Fr(1), Fr(7, 3), 4, 1), (Fr(2, 3), Fr(1), 1, 4), (Fr(5), Fr(2), 4, 2),
           (Fr(1), Fr(1, 9), 3, 3), (Fr(4, 3), Fr(1), 5, 1), (Fr(1), Fr(1), 4, 0), (Fr(2), Fr(1), 5, 0),
           (Fr(1), Fr(1, 2), 4, 3), (Fr(3, 2), Fr(1), 5, 2)]
    for (w0, w1, u0, u1) in pts:
        good_c, tot = count_check(u0, u1)
        ok &= good_c
        print(f'(u0,u1) = ({u0},{u1}): dimension count {tot} vs (k^2, k(k-1)/2): {good_c}', flush=True)
        for b in (Fr(1), Fr(1001, 1000)):
            for part in ('inv', 'std', 'two', 'anti', 'mix'):
                if part == 'std' and u0 < 2:
                    continue
                if part == 'two' and u0 < 4:
                    continue
                if part == 'anti' and u0 < 3:
                    continue
                if part == 'mix' and (u0 < 2 or u1 < 2):
                    continue
                kp, kq, gauge = keep_for(part, u0, u1)
                Mo = Model(Fr(1), [w0 * w0, w1 * w1], [Fr(u0), Fr(u1)], b)
                red = reduce_part(Mo, part, kp, kq)
                ex = explicit(part, w0, w1, u0, u1, b)
                same = red['RPP'] == ex['RPP'] and red['S'] == ex['S']
                if 'RPQ' in ex:
                    same &= red['RPQ'] == ex['RPQ'] and red['RQQ'] == ex['RQQ']
                gram_ok = ex['GP'] == red['GP']
                nonsing = det(ex['GP'], Fr(1)) != 0
                kerok = True
                if b == 1:
                    for kv in red['ker']:
                        kerok &= all(sum(red['S'][i][j] * kv[j] for j in range(len(kv))) == 0 for i in range(len(kv)))
                good = same and gram_ok and nonsing and kerok
                ok &= good
                print(f'  b = {b}, {part:4s} (P {len(kp)}, Q {len(kq) if kq else 0}): reduced = explicit {same}, '
                      f'Gram {gram_ok}, basis {nonsing}, kernel {kerok} [{time.time() - t0:.1f}s]', flush=True)
    print(f'time {time.time() - t0:.1f}s')
    print('PASS' if ok else 'FAIL')


if __name__ == '__main__':
    main()
