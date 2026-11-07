"""Measurements of the split (tsplit.py, tsbuild.py) on a sym family: the sizes of Q2g and Q3, the degree bounds of
T = D* det for each part (with the structure lemma when i and tau' are symbolic), the time of one slice, and the
symmetric elimination of Q3 (with a time limit).  Nothing here is part of a proof.
python tsmeasure.py FAM [--q3=elim|det|none] [--q2=det|none] [--tlimit=S] [--slice=0|1]   (FAM as in tprove.py;
FREE=size selects the free size of the slice when i and tau' are not both symbolic)"""
import os, sys, time
import flint
import tfamlib as T
import tprove as P
import tdet
import tsbuild


def to_nu_any(C, Qg, maxdeg):
    """as tdet.to_nu, allowing degree <= maxdeg in i"""
    old = C.syms
    new = [('nu' if nm == 'tp' else nm) for nm in old]
    ctx2 = flint.fmpq_mpoly_ctx.get(tuple(new), 'lex')
    g2 = dict(zip(new, ctx2.gens()))
    args = [g2['nu'] - 2 * g2['i'] if nm == 'tp' else g2[nm] for nm in old]
    ii = new.index('i')
    E, worst = {}, 0
    for a, row in enumerate(Qg):
        for b, x in enumerate(row):
            if x == 0:
                continue
            num = x.n.compose(*args, ctx=ctx2)
            den = x.d.compose(*args, ctx=ctx2)
            g = num.gcd(den)
            num, den = num / g, den / g
            assert den.degrees()[ii] == 0, 'denominator depends on i'
            worst = max(worst, num.degrees()[ii])
            E[(a, b)] = (num, den)
    assert worst <= maxdeg, f'degree {worst} in i'
    return E, ctx2, new, worst


def measure_part(C, Qg, label, log, maxdeg=2, slice_too=True):
    n = len(Qg)
    if n == 0:
        log(f'  {label}: empty')
        return
    T0 = time.time()
    structured = C.fixed['i'] is None and C.fixed['tp'] is None
    if structured:
        E, ctx2, names2, worst = to_nu_any(C, Qg, maxdeg)
        log(f'  {label}: structure holds (denominators free of i, degree {worst} in i at fixed nu) '
            f'[{time.time()-T0:.1f}s]')
    else:
        E = {(a, b): (x.n, x.d) for a, row in enumerate(Qg) for b, x in enumerate(row) if x != 0}
        ctx2, names2 = None, C.syms
    big = max(len(nm) for nm, dn in E.values())
    bigd = max(len(dn) for nm, dn in E.values())
    facs, estar = tdet.factor_data(E, n)
    bounds = {}
    for iv, v in enumerate(names2):
        if structured and v == 'i':
            W = [[(E[(a, b)][0].degrees()[iv] if (a, b) in E else None) for b in range(n)] for a in range(n)]
            bounds[v] = tdet.maxperm(W)
        else:
            bounds[v] = tdet.degree_bound(E, n, facs, estar, iv)
    log(f'  {label}: {n} x {n}, numerators up to {big} terms, denominators up to {bigd}; {len(facs)} irreducible '
        f'factors, exponents of D* {estar}; degree bounds of T = D* det: {bounds} [{time.time()-T0:.1f}s]')
    if not slice_too:
        return bounds
    if structured:
        fix = {v: C.smin[v] for v in names2 if v in ('q0', 't', 'i') and v in C.smin and C.fixed[v] is None}
        free = ['nu', 's']
    else:
        szs = [v for v in names2 if v != 's']
        fv = os.environ.get('FREE', szs[0])
        free = [fv, 's']
        fix = {v: C.smin[v] for v in szs if v != fv}
    t1 = time.time()
    K, det = tdet.slice_det(E, n, ctx2, names2, fix, free)
    log(f'  {label}: one slice ({fix} fixed, det over Q({", ".join(free)})): {len(det.n)}/{len(det.d)} terms, '
        f'{time.time()-t1:.1f}s')
    return bounds


if __name__ == '__main__':
    opts = {'q3': 'elim', 'q2': 'det', 'tlimit': '600', 'slice': '1'}
    args = []
    for a in sys.argv[1:]:
        if a.startswith('--'):
            k_, v_ = a[2:].split('=')
            opts[k_] = v_
        else:
            args.append(a)
    logf = lambda m: print(m, flush=True)
    for spec in args:
        case, mode, fam = P.parse_fam(spec)
        print(P.fam_label(case, mode, fam), flush=True)
        T0 = time.time()
        R = tsbuild.build_split(case, fam, logf)
        C = R['C']
        if opts['q3'] == 'elim' and R['three']:
            n3 = len(R['three'])
            M = {(a, b): R['Q3'][a][b] for a in range(n3) for b in range(n3)}
            res, info = P.eliminate(M, list(range(n3)), C, time.time(), float(opts['tlimit']), True, R['three'])
            print(f'  Q3 elimination: {res}; pivots {info["pivots"]} [{time.time()-T0:.1f}s]', flush=True)
        elif opts['q3'] == 'det' and R['three']:
            measure_part(C, R['Q3'], 'Q3', logf, maxdeg=4)
        if opts['q2'] == 'det':
            measure_part(C, R['Q2g'], 'Q2g', logf, maxdeg=2, slice_too=opts.get('slice', '1') == '1')
        print(f'  total {time.time()-T0:.1f}s', flush=True)
    print('PASS')
