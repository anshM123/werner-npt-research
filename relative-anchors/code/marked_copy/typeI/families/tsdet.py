"""Propositions D and D1 (NOTES.md, Sections 4 and 4b) on the first part of the split form, with the second part by
elimination (NOTES.md, Section 8).

The builder replaces tdet.build: tsbuild.build_split (preliminaries, the change of basis with its determinant one-signed,
the block-diagonality T^T Q T = Q2 (+) Q3 identically, the kernel vectors in the first part annihilated, the gauge on the
first part), then Q3 by greedy symmetric elimination with pivots one-signed after the shifts (tprove.eliminate); the
gauge-fixed first part Q2g is handed to tdet.prove (Proposition D) or tdet1.prove1 (Proposition D1) unchanged.  At every
member and every sigma > 0, Q >= 0 with kernel exactly span(kappa) iff Q2g > 0 and Q3 > 0.
python tsdet.py D OUT.pkl CASE A,C [A,C ...]      (the families [q0 = A, t = C, i >= 2, tau' >= 2])
python tsdet.py D1 OUT.pkl FAM [FAM ...]          (two symbolic sizes, not both i and tau')"""
import os, sys, time, pickle
import tfamlib as T
import tprove as P
import tdet
import tdet1
import tsbuild

ORDER2 = {'MT': 'U0,U1,W0,W1', 'MU': 'U00,U01,U10,U11,W00,W01,W10,W11'}


def build_split_for_det(case, fam, log=print):
    os.environ['ORDER2'] = os.environ.get('ORDER2_' + case, ORDER2[case])
    R = tsbuild.build_split(case, fam, log)
    C = R['C']
    n3 = len(R['three'])
    if n3:
        t0 = time.time()
        M = {(a, b): R['Q3'][a][b] for a in range(n3) for b in range(n3)}
        res, info = P.eliminate(M, list(range(n3)), C, t0, 3600, False, R['three'])
        assert res == 'OK', ('Q3 elimination', res)
        log(f'  Q3 ({n3} x {n3}) positive definite: pivots {info["pivots"]} (one-signed after the shifts) '
            f'[{time.time()-t0:.1f}s]')
    else:
        log('  Q3 empty (no split block)')
    return C, R['names'], R['keep2'], R['Q2g']


tdet.build = build_split_for_det


if __name__ == '__main__':
    mode, out = sys.argv[1], sys.argv[2]
    logf = lambda m: print(m, flush=True)
    T0 = time.time()
    results, infos = {}, {}
    allok = True
    if mode == 'D':
        case = sys.argv[3]
        for spec in sys.argv[4:]:
            q0, t = (int(x) for x in spec.split(','))
            fam = (q0, t, None, None)
            lab = P.fam_label(case, 'sym', fam)
            print(lab, flush=True)
            key = (case, 'sym', q0, t, ('>=', T.SMIN['i']), ('>=', T.SMIN['tp']))
            try:
                ok, info = tdet.prove(case, q0, t, None, logf)
                res = 'PD on complement (determinant route with the split)' if ok else 'FAILED (a sign or the base point)'
            except AssertionError as ex:
                ok, info, res = False, {'time': 0}, f'FAILED ({ex})'
            allok &= ok
            results[key], infos[key] = res, info
            print(f'{lab} {res} ({info["time"]:.1f}s)', flush=True)
    else:
        for s in sys.argv[3:]:
            case, mode2, fam = P.parse_fam(s)
            assert mode2 == 'sym' and sum(x is None for x in fam) == 2 and not (fam[2] is None and fam[3] is None)
            lab = P.fam_label(case, mode2, fam)
            print(lab, flush=True)
            key = (case, mode2) + tuple(('>=', T.SMIN[nm]) if x is None else x for x, nm in zip(fam, T.SIZES))
            try:
                ok, info = tdet1.prove1(case, fam, logf)
                res = ('PD on complement (determinant route with the split, one size sliced)' if ok
                       else 'FAILED (a sign or the base point)')
            except AssertionError as ex:
                ok, info, res = False, {'time': 0}, f'FAILED ({ex})'
            allok &= ok
            results[key], infos[key] = res, info
            print(f'{lab} {res} ({info["time"]:.1f}s)', flush=True)
    with open(out, 'wb') as fh:
        pickle.dump({'results': results, 'info': infos, 'smin': dict(T.SMIN)}, fh)
    print(f'{sum(r.startswith("PD") for r in results.values())} families proved of {len(results)}; '
          f'{time.time() - T0:.0f}s')
    print('PASS' if allok else 'FAIL')
