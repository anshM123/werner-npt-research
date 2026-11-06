"""Family certificates with the irrep split (splitlib.py).  For a family (q0, t, i each fixed or symbolic; symbolic
means q0 = 4 + x, t = 4 + y, i = 2 + z with x, y, z >= 0; sigma > 0 symbolic):
  (1) the reduced form Q over Q(sizes, sigma) (famlib.reduced), and the sanity check of prove_fam (N0 times Q equals
      the exact numeric form of mlib at 3^s sample points and two values of sigma; mlib asserts its Newton window);
  (2) every entry denominator is one-signed after the shifts (so Q is the true form at every member and sigma > 0);
  (3) the change of basis T of splitlib (on each split block det = m + 1 + i, one-signed after the shifts) and the
      identity T^T Q T = Q2 (+) Q3 (the off-diagonal block vanishes identically);
  (4) the known kernel vectors of prove_fam.fam_kernel, written in the new basis (they lie in the first part), are
      annihilated identically; gauge: delete one coordinate per kernel vector, with a one-signed determinant;
  (5) greedy symmetric elimination of each part with pivots that are ratios of one-signed polynomials after the
      shifts (as in prove_fam: then every leading principal minor of a symmetric permutation is > 0 on the family).
Then Q >= 0 with kernel exactly the span of the known kernel vectors, at every member and every sigma > 0.
Usage: python prove_split.py q0 t i [bscale] [tlimit]   ('n' for a symbolic size)"""
import sys, time, itertools, pickle
from fractions import Fraction as Fr
import famlib, mlib, splitlib
from prove_fam import one_signed, evalF, SMIN, fam_kernel, det_small
import os
GAUGE = os.environ.get('GAUGE', 'O')   # 'O': delete O_{y0}; 'R': delete U_{1,y0} (A_{1,y0})


def sanity(C, names, Q, q0, t, i, bscale):
    sz = [nm for nm in C.syms if nm != 's']
    for sample in itertools.product(*[(SMIN[nm], SMIN[nm] + 1, SMIN[nm] + 3) for nm in sz]):
        val = dict(zip(sz, sample))
        q0v, tv, iv = val.get('q0', q0), val.get('t', t), val.get('i', i)
        m = q0v + tv
        u = m + 2 * iv
        k, j = u + 1, q0v + iv
        N0 = famlib.norm0(C, q0v, tv, iv)
        for sg in (Fr(2, 7), Fr(5)):
            n2, nr2, Q2 = mlib.reduced(mlib.Calc(k, j, iv, sigma=sg, bscale=bscale))
            if n2 != names:
                return f'SANITY-NAMES {names} {n2} at {(k, j, iv)}'
            pt = list(sample) + [sg]
            for a in range(len(names)):
                for b in range(len(names)):
                    if evalF(Q[a][b], pt) * N0 != Fr(int(Q2[a][b].p), int(Q2[a][b].q)):
                        return f'SANITY-VALUES at {(k, j, iv, sg)} {names[a]},{names[b]}'
    return None


def eliminate(M, act, C, label, verbose, t0, tlimit=None):
    """greedy symmetric elimination on the name list act; returns (status, info)"""
    act = list(act)
    info = {'order': [], 'maxterms': 0, 'pivots': []}
    while act:
        if tlimit is not None and time.time() - t0 > tlimit:
            return 'TIME', info
        cands = sorted(act, key=lambda z: (len(M[z, z].n) + len(M[z, z].d)) if M[z, z] != 0 else 10 ** 12)
        chosen = None
        for x in cands:
            piv = M[x, x]
            if piv == 0:
                continue
            sn, sd = one_signed(piv.n, C), one_signed(piv.d, C)
            if sn != 0 and sd != 0 and sn == sd:
                chosen = x
                break
        if chosen is None:
            x = cands[0]
            return (f'PIVOT-SIGN-UNPROVED in {label} (smallest {x}: num {len(M[x, x].n)} terms, '
                    f'den {len(M[x, x].d)})'), info
        x = chosen
        piv = M[x, x]
        info['maxterms'] = max(info['maxterms'], len(piv.n), len(piv.d))
        info['order'].append(x)
        info['pivots'].append((x, len(piv.n), len(piv.d)))
        if verbose:
            print(f'    {label} pivot {x}: num {len(piv.n)} terms, den {len(piv.d)} terms ({time.time()-t0:.1f}s)',
                  flush=True)
        act.remove(x)
        for ia, a in enumerate(act):
            if M[a, x] == 0:
                continue
            f = M[a, x] / piv
            for b in act[ia:]:
                if M[x, b] != 0:
                    M[a, b] = M[a, b] - f * M[x, b]
                    M[b, a] = M[a, b]
    return 'OK', info


def prove_split(q0, t, i, verbose=False, bscale=Fr(1), tlimit=None):
    t0 = time.time()
    if q0 is not None and t is not None and q0 + t < 0:
        return 'empty', 0, {}
    C = famlib.FamCalc(q0, t, i, bscale=bscale)
    names, Q = famlib.reduced(C)
    if not names:
        return 'empty', time.time() - t0, {}
    if verbose:
        print(f'  form computed ({time.time()-t0:.1f}s), {len(names)} vectors', flush=True)
    bad = sanity(C, names, Q, q0, t, i, bscale)
    if bad:
        return bad, time.time() - t0, {}
    n = len(names)
    for a in range(n):
        for b in range(n):
            if Q[a][b] != 0 and one_signed(Q[a][b].d, C) == 0:
                return f'ENTRY DENOMINATOR NOT ONE-SIGNED at {names[a]},{names[b]}', time.time() - t0, {}
    vecs, split = splitlib.basis(names, C.q0, C.m, C.i, C.one)
    if split:
        d = C.m + 1 + C.i
        if one_signed(d.n, C) != 1 or one_signed(d.d, C) != 1:
            return 'DET T NOT POSITIVE', time.time() - t0, {}
    two, three, M = splitlib.congruence(names, Q, vecs, C.zero)
    if verbose:
        print(f'  congruence done ({time.time()-t0:.1f}s); split blocks {split}; parts {len(two)} + {len(three)}',
              flush=True)
    for a in two:
        for b in three:
            if M[a, b] != 0:
                return f'SPLIT NOT BLOCK DIAGONAL at {a},{b}', time.time() - t0, {}
    kv = []
    fk = fam_kernel(names, C)
    # gauge R only when kappa_0, kappa_1 are the only known kernel vectors (else the coordinates may clash)
    gauge_r = GAUGE == 'R' and all(lab.startswith('p x h') for lab, _, _ in fk)
    for lab, v, dl in fk:
        w = splitlib.to_new(v, split)
        if gauge_r and w is not None:
            # delete the row-in-U coordinate U_{1,y0} (or A_{1,y0} if that block is not split): coefficient 1
            y0 = lab[-1]
            for cand in ('U1' + y0, 'A1' + y0):
                if cand in w and w[cand] == 1:
                    dl = cand
                    break
        if w is not None and dl not in w:
            # the designated coordinate was split away: take the first coordinate of the first part, not yet used,
            # where the vector has a nonzero constant coefficient (the gauge determinant is checked below)
            used = {d for _, _, d in kv}
            const = lambda c: c != 0 and all(sum(e) == 0 for e in c.n.to_dict()) and                 all(sum(e) == 0 for e in c.d.to_dict())
            cands = [nm for nm in two if nm in w and nm not in used and const(w[nm])]
            if cands:
                dl = cands[0]
        if w is None or dl not in w or any(nm not in two for nm in w):
            return f'KERNEL VECTOR {lab} NOT IN THE FIRST PART OR GAUGE COORDINATE SPLIT', time.time() - t0, {}
        kv.append((lab, w, dl))
    if bscale == 1:
        for lab, v, _ in kv:
            for a in two:
                s = C.zero
                for nm, cf in v.items():
                    s = s + M[a, nm] * cf
                if s != 0:
                    return f'KAPPA {lab} NOT IN KERNEL', time.time() - t0, {}
    gauge = [dl for _, _, dl in kv]
    sub = [[v.get(c, C.zero) for c in gauge] for _, v, _ in kv]
    det = det_small(sub, C)
    if det == 0 or one_signed(det.n, C) == 0 or one_signed(det.d, C) == 0:
        return f'GAUGE DETERMINANT NOT SIGN-DEFINITE', time.time() - t0, {}
    info = {'names': names, 'two': two, 'three': three, 'gauge': gauge, 'kernel': [lab for lab, _, _ in kv],
            'split': split}
    res3, info3 = eliminate(M, three, C, 'part3', verbose, t0, tlimit)
    info['part3'] = info3
    if res3 != 'OK':
        return res3, time.time() - t0, info
    keep2 = [nm for nm in two if nm not in gauge]
    res2, info2 = eliminate(M, keep2, C, 'part2', verbose, t0, tlimit)
    info['part2'] = info2
    if res2 != 'OK':
        return res2, time.time() - t0, info
    info['maxterms'] = max(info2['maxterms'], info3['maxterms'])
    return (f'PD on complement (pivots {len(keep2)} + {len(three)}, kernel {len(kv)})', time.time() - t0, info)


def parse(s):
    return None if s == 'n' else int(s)


if __name__ == '__main__':
    q0, t, i = (parse(s) for s in sys.argv[1:4])
    bscale = Fr(sys.argv[4]) if len(sys.argv) > 4 else Fr(1)
    tl = float(sys.argv[5]) if len(sys.argv) > 5 else None
    res, dt, info = prove_split(q0, t, i, verbose=True, bscale=bscale, tlimit=tl)
    print(f'[{q0}, {t}, {i}] {res} ({dt:.1f}s)', flush=True)
    print('PASS' if res.startswith('PD') else 'FAIL')
