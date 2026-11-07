"""The split family form (tsplit.py) with the preliminaries of tprove.py / tdet.build.

build_split(case, fam) for a sym family (each of q0, t, i, tau' fixed or symbolic):
  1. the form Q over Q(sizes, sigma) (tfamlib.reduced), symmetric; N0 Q equals tlib's numeric form at the sample points
     (each symbolic size at min, min + 1, min + 3; sigma = 2/7, 5); every intermediate denominator of the calculus and
     every entry denominator one-signed after the shifts (Lemma V: Q is the true form at every member, sigma > 0);
  2. the change of basis of tsplit.py: on each split block its determinant m + 1 + i is one-signed after the shifts,
     and T^T Q T = Q2 (+) Q3 identically (the block between the two parts vanishes as rational functions);
  3. the known kernel vectors (tlib.kernel_vectors), written in the new basis, have no component on the S vectors and
     are annihilated by Q2 identically;
  4. gauge on Q2: one deleted coordinate per kernel vector (preferred order GPREF), the determinant of the kernel vectors
     there a nonzero constant.
Then at every member and every sigma > 0: Q >= 0 with kernel exactly span(kappa) iff the gauge-fixed Q2g is positive
definite and Q3 is positive definite (T is invertible, the kernel lies in the first part, and the gauge coordinates
give a complement of it there).
Returns a dict with C, the old names, the two-part names kept (keep2), the gauge, Q2g (list of lists), the
three-part names and Q3 (list of lists)."""
import os, time
from fractions import Fraction as Fr
import tfamlib as T
import tprove as P
import tsplit

GPREF = {'MT': os.environ.get('GAUGE_MT', 'CM,C,U1').split(','),
         'MU': os.environ.get('GAUGE_MU', 'O0,O1,C0,C1,U10,U11').split(',')}


def build_split(case, fam, log=print, do_sanity=True, bscale=1):
    t0 = time.time()
    C = T.FamCalc(case, 'sym', *fam, bscale=bscale)
    assert C.valid(), 'family not valid'
    names, Q = T.reduced(C)
    n = len(names)
    assert n > 0, 'empty family'
    log(f'  form: {n} vectors, largest entry {max(P.tsz(Q[a][b]) for a in range(n) for b in range(n))} terms '
        f'[{time.time()-t0:.1f}s]')
    assert all(Q[a][b] == Q[b][a] for a in range(n) for b in range(n)), 'not symmetric'
    npts = 0
    if do_sanity:
        bad, npts = P.sanity(C, names, Q, Fr(bscale))
        assert bad is None, bad
    for key, (x, why) in C.dens.items():
        assert P.one_signed(x.n, C) != 0 and x.d.is_constant(), ('intermediate denominator', why)
    for a in range(n):
        for b in range(n):
            assert Q[a][b] == 0 or P.one_signed(Q[a][b].d, C) != 0, ('entry denominator', names[a], names[b])
    vecs, two, three, split = tsplit.basis(case, names, C.q0, C.m, C.i, C.one)
    if split:
        dT = C.m + 1 + C.i
        assert P.one_signed(dT.n, C) == 1 and dT.d.is_constant(), 'change of basis determinant not one-signed'
    order = two + three
    M = tsplit.congruence(names, Q, vecs, order, C.zero)
    nz = [(p, q) for p in two for q in three if M[p, q] != 0]
    assert not nz, ('the split is not block-diagonal', nz[:3])
    log(f'  split: {len(split)} split blocks, parts {len(two)} + {len(three)}, block-diagonal identically '
        f'[{time.time()-t0:.1f}s]')
    kv = T.kernel_vectors(C, names)
    kv2 = []
    for lab, v in kv:
        w = tsplit.to_new(v, split)
        assert w is not None and not any(nm in three for nm in w), ('kernel vector not in the first part', lab)
        for p in two:
            s = C.zero
            for nm, cf in w.items():
                s = s + M[p, nm] * cf
            assert s == 0, ('kernel vector not annihilated', lab, p)
        kv2.append((lab, w))
    # gauge on the first part
    r = len(kv2)
    gauge = []
    if r:
        import itertools
        pref = [nm for nm in GPREF[case] if nm in two] + [nm for nm in two if nm not in GPREF[case]]
        for S in itertools.combinations(pref, r):
            sub = [[C.K(v.get(nm, 0)) for nm in S] for _, v in kv2]
            d = P.det_small(sub, C)
            if d != 0 and d.n.is_constant() and d.d.is_constant():
                gauge = list(S)
                break
        assert gauge, 'no gauge'
    keep2 = [nm for nm in two if nm not in gauge]
    first = [nm for nm in os.environ.get('ORDER2', '').split(',') if nm in keep2]
    keep2 = first + [nm for nm in keep2 if nm not in first]
    Q2g = [[M[a, b] for b in keep2] for a in keep2]
    Q3 = [[M[a, b] for b in three] for a in three]
    big2 = max((P.tsz(x) for row in Q2g for x in row), default=0)
    big3 = max((P.tsz(x) for row in Q3 for x in row), default=0)
    log(f'  preliminaries: symmetric, {npts} sanity points, {len(C.dens)} intermediate denominators and the entry '
        f'denominators one-signed, {r} kernel vectors in the first part annihilated, gauge {gauge}; Q2g {len(keep2)} x '
        f'{len(keep2)} (largest entry {big2} terms), Q3 {len(three)} x {len(three)} (largest entry {big3} terms) '
        f'[{time.time()-t0:.1f}s]')
    return {'C': C, 'names': names, 'two': two, 'three': three, 'keep2': keep2, 'gauge': gauge, 'Q2g': Q2g,
            'Q3': Q3, 'kv': kv2, 'split': split, 'npts': npts}
