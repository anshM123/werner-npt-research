"""Development probe (no proof): time of one slice determinant of the split first part for several orders of the
elimination (tdet.slice_det eliminates in the order of the rows)."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import tprove as P
import tdet, tsbuild

spec, fixs = sys.argv[1], sys.argv[2]
case, mode, fam = P.parse_fam(spec)
orders = sys.argv[3:]
fix = {kv.split('=')[0]: int(kv.split('=')[1]) for kv in fixs.split(';')}
for od in orders:
    os.environ['ORDER2'] = od
    R = tsbuild.build_split(case, fam, lambda m: None, do_sanity=False)
    C = R['C']
    structured = fam[2] is None and fam[3] is None
    if structured:
        E, ctx2, names2 = tdet.to_nu(C, R['Q2g'])
        free = ['nu', 's']
    else:
        E = {(a, b): (x.n, x.d) for a, row in enumerate(R['Q2g']) for b, x in enumerate(row) if x != 0}
        ctx2, names2 = None, C.syms
        free = [v for v in names2 if v not in fix]
    t1 = time.time()
    K, det = tdet.slice_det(E, len(R['keep2']), ctx2, names2, fix, free)
    print(f'{od or "(default)"}: keep {R["keep2"]}: det {len(det.n)}/{len(det.d)} terms, {time.time()-t1:.2f}s', flush=True)


def slice_det_dyn(E, n, ctx2, names2, fix, free):
    """the determinant with symmetric pivoting on the smallest remaining diagonal entry (number of terms)"""
    from flintfield import Field, F
    K = Field(list(free))
    gens = dict(zip(free, K.ctx.gens()))
    args = [gens[v] if v in gens else K.ctx.from_dict({(0,) * len(free): fix[v]}) for v in names2]
    A = {(a, b): F(K, E[(a, b)][0].compose(*args, ctx=K.ctx), E[(a, b)][1].compose(*args, ctx=K.ctx))
         for (a, b) in E}
    act = list(range(n))
    det = K(1)
    while act:
        cand = [x for x in act if A.get((x, x), 0) != 0]
        x = min(cand, key=lambda z: len(A[(z, z)].n) + len(A[(z, z)].d))
        p = A[(x, x)]
        det = det * p
        act.remove(x)
        for r in act:
            ar = A.get((r, x), 0)
            if ar == 0:
                continue
            f = ar / p
            for c in act:
                xc = A.get((x, c), 0)
                if xc != 0:
                    A[(r, c)] = A.get((r, c), K(0)) - f * xc
    return K, det


if os.environ.get('DYN') == '1':
    os.environ['ORDER2'] = ''
    R = tsbuild.build_split(case, fam, lambda m: None, do_sanity=False)
    C = R['C']
    structured = fam[2] is None and fam[3] is None
    if structured:
        E, ctx2, names2 = tdet.to_nu(C, R['Q2g'])
        free = ['nu', 's']
    else:
        E = {(a, b): (x.n, x.d) for a, row in enumerate(R['Q2g']) for b, x in enumerate(row) if x != 0}
        ctx2, names2 = None, C.syms
        free = [v for v in names2 if v not in fix]
    t1 = time.time()
    K, det = slice_det_dyn(E, len(R['keep2']), ctx2, names2, fix, free)
    print(f'dynamic smallest pivot: det {len(det.n)}/{len(det.d)} terms, {time.time()-t1:.2f}s', flush=True)
