"""Family certificates with the irrep split (as prove_split.py), for families whose symbolic sizes have chosen lower
bounds: q0 = QMIN + x, t = TMIN + y, i = IMIN + z (defaults 4, 4, 5), sigma > 0 symbolic.  The bounds are written
into the shared dictionary famlib.SMIN before anything is computed, so the presence of the twisted vectors, the
sanity points (each symbolic size at its minimum, minimum + 1, minimum + 3) and the shifts of the sign test all use
them.  The formal continuation of famlib.py is valid for symbolic q0 >= 4, t >= 4, i >= 2, which these bounds imply.
Steps, as in prove_split.py: (1) the form over Q(sizes, sigma) and the sanity check against mlib.py (which asserts
its Newton window); (2) every entry denominator one-signed after the shifts; (3) the change of basis of splitlib.py
(det m + 1 + i per split block, one-signed) and T^T Q T = Q2 (+) Q3 identically; (4) the known kernel vectors in the
first part, annihilated identically; gauge deleting U_{1,y0} (or A_{1,y0}), determinant one-signed; (5) greedy
symmetric elimination of the parts asked for (--parts=3: the (u - i - 1, i, 1) part only; --parts=2: the
(u - i, i) part only; default both), every pivot a ratio of one-signed polynomials after the shifts.
With --parts=3 the conclusion is only that Q3 (the form on the S vectors) is positive definite.
Usage: python prove_gen.py OUT.pkl FAM [FAM ...] [--imin=5] [--qmin=4] [--tmin=4] [--parts=32] [--bscale=1]
       [--tlimit=SECONDS]      (FAM = q0,t,i with 'n' for a symbolic size; OUT.pkl is written only when bscale = 1)
The saved keys are (q0, t, i) with a symbolic size written ('>=', minimum)."""
import sys, os, time, pickle
from fractions import Fraction as Fr

opts = {'imin': '5', 'qmin': '4', 'tmin': '4', 'parts': '32', 'bscale': '1', 'tlimit': ''}
args = []
for a in sys.argv[1:]:
    if a.startswith('--'):
        k, v = a[2:].split('=')
        opts[k] = v
    else:
        args.append(a)

import famlib
famlib.SMIN['q0'] = int(opts['qmin'])
famlib.SMIN['t'] = int(opts['tmin'])
famlib.SMIN['i'] = int(opts['imin'])
assert famlib.SMIN['q0'] >= 4 and famlib.SMIN['t'] >= 4 and famlib.SMIN['i'] >= 2
import splitlib
from prove_fam import one_signed, fam_kernel, det_small, SMIN
from prove_split import sanity
assert SMIN is famlib.SMIN


def tsz(x):
    return 0 if x == 0 else len(x.n) + len(x.d)


def eliminate(M, act, C, label, t0, tlimit=None):
    """greedy symmetric elimination (as prove_split.eliminate) with progress lines"""
    act = list(act)
    info = {'order': [], 'maxterms': 0, 'pivots': []}
    while act:
        if tlimit is not None and time.time() - t0 > tlimit:
            return 'TIME', info
        cands = sorted(act, key=lambda z: tsz(M[z, z]) if M[z, z] != 0 else 10 ** 12)
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
        info['pivots'].append((x, len(piv.n), len(piv.d), round(time.time() - t0, 1)))
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
            if len(act) > 1:
                big = max(tsz(M[a, b]) for b in act)
                print(f'      row {a} updated ({time.time()-t0:.1f}s), largest entry {big} terms', flush=True)
        if tlimit is not None and time.time() - t0 > tlimit:
            return 'TIME', info
    return 'OK', info


def prove_gen(q0, t, i, parts='32', bscale=Fr(1), tlimit=None):
    t0 = time.time()
    if q0 is not None and t is not None and q0 + t < 0:
        return 'empty', 0, {}
    C = famlib.FamCalc(q0, t, i, bscale=bscale)
    names, Q = famlib.reduced(C)
    if not names:
        return 'empty', time.time() - t0, {}
    big = max(tsz(Q[a][b]) for a in range(len(names)) for b in range(len(names)))
    print(f'  form computed ({time.time()-t0:.1f}s), {len(names)} vectors, largest entry {big} terms', flush=True)
    bad = sanity(C, names, Q, q0, t, i, bscale)
    if bad:
        return bad, time.time() - t0, {}
    print(f'  sanity against mlib done ({time.time()-t0:.1f}s)', flush=True)
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
    big2 = max([tsz(M[a, b]) for a in two for b in two] or [0])
    big3 = max([tsz(M[a, b]) for a in three for b in three] or [0])
    print(f'  congruence done ({time.time()-t0:.1f}s); split blocks {split}; parts {len(two)} + {len(three)}; '
          f'largest entries {big2} / {big3} terms', flush=True)
    for a in two:
        for b in three:
            if M[a, b] != 0:
                return f'SPLIT NOT BLOCK DIAGONAL at {a},{b}', time.time() - t0, {}
    kv = []
    fk = fam_kernel(names, C)
    gauge_r = all(lab.startswith('p x h') for lab, _, _ in fk)
    for lab, v, dl in fk:
        w = splitlib.to_new(v, split)
        if gauge_r and w is not None:
            y0 = lab[-1]
            for cand in ('U1' + y0, 'A1' + y0):
                if cand in w and w[cand] == 1:
                    dl = cand
                    break
        if w is not None and dl not in w:
            used = {d for _, _, d in kv}
            const = lambda c: c != 0 and all(sum(e) == 0 for e in c.n.to_dict()) and \
                all(sum(e) == 0 for e in c.d.to_dict())
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
        return 'GAUGE DETERMINANT NOT SIGN-DEFINITE', time.time() - t0, {}
    info = {'names': names, 'two': two, 'three': three, 'gauge': gauge, 'kernel': [lab for lab, _, _ in kv],
            'split': split, 'parts': parts}
    info['maxterms'] = 0
    if '3' in parts:
        res3, info3 = eliminate(M, three, C, 'part3', t0, tlimit)
        info['part3'] = info3
        info['maxterms'] = max(info['maxterms'], info3['maxterms'])
        if res3 != 'OK':
            return res3, time.time() - t0, info
    keep2 = [nm for nm in two if nm not in gauge]
    if '2' in parts:
        res2, info2 = eliminate(M, keep2, C, 'part2', t0, tlimit)
        info['part2'] = info2
        info['maxterms'] = max(info['maxterms'], info2['maxterms'])
        if res2 != 'OK':
            return res2, time.time() - t0, info
    if parts == '3':
        return f'PART3 PD (pivots {len(three)})', time.time() - t0, info
    if parts == '2':
        return f'PART2 PD on complement (pivots {len(keep2)}, kernel {len(kv)})', time.time() - t0, info
    return (f'PD on complement (pivots {len(keep2)} + {len(three)}, kernel {len(kv)})', time.time() - t0, info)


def parse(s):
    return None if s == 'n' else int(s)


def key_of(fam):
    return tuple(('>=', SMIN[nm]) if v is None else v for v, nm in zip(fam, ('q0', 't', 'i')))


def label(fam):
    return '[' + ', '.join(f'{nm}>={SMIN[nm]}' if v is None else f'{nm}={v}'
                           for v, nm in zip(fam, ('q0', 't', 'i'))) + ']'


if __name__ == '__main__':
    out = args[0]
    fams = [tuple(parse(s) for s in a.split(',')) for a in args[1:]]
    bscale = Fr(opts['bscale'])
    tlimit = float(opts['tlimit']) if opts['tlimit'] else None
    parts = opts['parts']
    T0 = time.time()
    allok = True
    results = {}
    infos = {}
    nproved = nempty = 0
    for fam in fams:
        print(f'{label(fam)} parts {parts}, bscale {bscale}', flush=True)
        res, dt, info = prove_gen(*fam, parts=parts, bscale=bscale, tlimit=tlimit)
        results[key_of(fam)] = res
        infos[key_of(fam)] = {k: v for k, v in info.items() if k in ('part2', 'part3', 'gauge', 'split', 'two',
                                                                     'three', 'names', 'maxterms')}
        ok = res.startswith('PD') or res.startswith('PART') or res == 'empty'
        allok &= ok
        nproved += res.startswith('PD') or res.startswith('PART')
        nempty += res == 'empty'
        print(f'{label(fam)} {res}  ({dt:.1f}s)' + (f'  largest pivot {info["maxterms"]} terms'
                                                    if info.get('maxterms') else ''), flush=True)
    if bscale == 1:
        with open(out, 'wb') as fh:
            pickle.dump({'results': results, 'parts': parts, 'info': infos}, fh)
    print(f'{nproved} families proved ({"parts " + parts}), {nempty} empty, of {len(fams)}; {time.time() - T0:.0f}s')
    print('PASS' if allok else 'FAIL')
