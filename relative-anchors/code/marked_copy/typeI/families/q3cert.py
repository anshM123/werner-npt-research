"""The three-row parts (NOTES.md, Section 8.6): for each sym family, the preliminaries and the split of tsbuild.py
(sanity against tlib.py at the sample points, every denominator one-signed after the shifts, the change of basis with
determinant m + 1 + i one-signed, T^T Q T = Q2 (+) Q3 identically, the kernel vectors in the first part), then Q3 by
greedy symmetric elimination with every pivot a ratio of one-signed polynomials of the same sign (tprove.eliminate).
Then at every member and every sigma > 0, Q3 is positive definite: the block is positive definite on the isotypic part
of the S_{U'}-irrep (u' - i - 1, i, 1) (for the members with a split block; Q3 is empty otherwise).
python q3cert.py OUT.pkl FAM [FAM ...]     (FAM as in tprove.py; @LIST reads a list)"""
import sys, time, pickle
import tfamlib as T
import tprove as P
import tsbuild


def cert(case, fam, log=print, tlimit=3600):
    t0 = time.time()
    R = tsbuild.build_split(case, fam, log)
    C = R['C']
    n3 = len(R['three'])
    if not n3:
        return 'no split block (Q3 empty)', time.time() - t0, {'three': []}
    M = {(a, b): R['Q3'][a][b] for a in range(n3) for b in range(n3)}
    res, info = P.eliminate(M, list(range(n3)), C, time.time(), tlimit, False, R['three'])
    info['three'] = R['three']
    if res != 'OK':
        return f'Q3 NOT CERTIFIED ({res})', time.time() - t0, info
    return f'Q3 positive definite ({n3} pivots)', time.time() - t0, info


if __name__ == '__main__':
    out = sys.argv[1]
    fams = []
    for a in sys.argv[2:]:
        if a.startswith('@'):
            fams += [P.parse_fam(s) for s in open(a[1:]).read().split()]
        else:
            fams.append(P.parse_fam(a))
    T0 = time.time()
    results, infos = {}, {}
    allok = True
    logf = lambda m: print(m, flush=True)
    for case, mode, fam in fams:
        lab = P.fam_label(case, mode, fam)
        print(lab, flush=True)
        try:
            res, dt, info = cert(case, fam, logf)
        except AssertionError as ex:
            res, dt, info = f'FAILED ({ex})', 0, {}
        key = (case, mode) + tuple(('>=', T.SMIN[nm]) if x is None else x for x, nm in zip(fam, T.SIZES))
        results[key], infos[key] = res, info
        ok = res.startswith('Q3 positive') or res.startswith('no split')
        allok &= ok
        piv = [(nm, a, b) for nm, a, b, _ in info.get('pivots', [])]
        print(f'{lab} {res}  ({dt:.1f}s) pivots (name, numerator terms, denominator terms): {piv}', flush=True)
    with open(out, 'wb') as fh:
        pickle.dump({'results3': results, 'info': infos}, fh)
    print(f'{sum(r.startswith("Q3 positive") for r in results.values())} families with Q3 certified, '
          f'{sum(r.startswith("no split") for r in results.values())} without a split block, of {len(results)}; '
          f'{time.time() - T0:.0f}s')
    print('PASS' if allok else 'FAIL')
