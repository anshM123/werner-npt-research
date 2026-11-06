"""The determinant route with slices in one size, for sym families with two symbolic sizes u, v (not both i and tau'),
as relative-anchors/code/middle_j2/f2det.py.  Qg is the gauge-fixed form (tdet.build: the preliminaries of tprove.py).  Steps:
 1. D* = prod f^{e_f} over the irreducible factors f of the entry denominators (e_f by dynamic programming over
    subsets); T = D* det Qg is a polynomial with deg_v T <= B = deg_v D* + max_p sum_a (deg_v num - deg_v den);
 2. slices v = vmin + y, y = 0..B + 1: det Qg exactly over Q(u, sigma); D* det must be a polynomial P_y on every slice
    (exact division);
 3. interpolation in y (inverse Vandermonde over Q): N~(y, u, sigma) of degree <= B with N~(y) = P_y for y = 0..B
    (through the extra node the coefficient of y^(B+1) must vanish, a consistency check of B).  T - N~ has degree <= B
    in v and vanishes at B + 1 values, so det Qg = N~/D* identically;
 4. signs: after the shifts every coefficient of N~ is >= 0 and one free of the size variables is > 0, and every factor
    of D* is one-signed positive after the shifts: det Qg > 0 on the region;
 5. Qg positive definite at the minimal member, sigma = 1 (exact elimination).
The region is connected and Qg is continuous there (one-signed entry denominators), so Qg > 0 on the whole region.
python tdet1.py OUT.pkl FAM [FAM ...]     (FAM as in tprove.py, sym mode, exactly two symbolic sizes)"""
import sys, time, pickle
import flint
from flintfield import Field
import tfamlib as T
import tprove as P
import tdet


def prove1(case, fam, log=print):
    T0 = time.time()
    C, names, keep, Qg = tdet.build(case, fam, log)
    n = len(keep)
    E = {(a, b): (x.n, x.d) for a, row in enumerate(Qg) for b, x in enumerate(row) if x != 0}
    names2 = C.syms
    szs = [v for v in names2 if v != 's']
    assert len(szs) == 2, szs
    facs, estar = tdet.factor_data(E, n)
    bounds = {v: tdet.degree_bound(E, n, facs, estar, names2.index(v)) for v in szs}
    v = min(szs, key=lambda w: bounds[w])          # slice in the size with the smaller bound
    u = [w for w in szs if w != v][0]
    iv = names2.index(v)
    B = bounds[v]
    log(f'  {len(facs)} irreducible factors, exponents of D* {estar}; bounds {bounds}; slices in {v}: 0..{B} '
        f'[{time.time()-T0:.1f}s]')
    free = [u, 's']
    vmin = C.smin[v]
    K = Field(free)
    gens = dict(zip(free, K.ctx.gens()))

    def spec(f, c):
        return f.compose(*[gens[w] if w in gens else K.ctx.from_dict({(0, 0): c}) for w in names2], ctx=K.ctx)

    # D_c = D* (all exponents e_f): R = 1 and T = D* det Qg is interpolated directly; B + 2 slices, the last one only
    # for a consistency check (the coefficient of y^(B+1) of the interpolant through B + 2 nodes must vanish)
    slices = []
    for y in range(B + 2):
        t1 = time.time()
        K2, det = tdet.slice_det(E, n, None, names2, {v: vmin + y}, free)
        assert det != 0, 'zero determinant on a slice'
        slices.append(det)
        if y % 20 == 0:
            log(f'    slice {v} = {vmin + y}: det {len(det.n)}/{len(det.d)} terms [{time.time()-t1:.1f}s]')
    obs = list(estar)
    L = B
    Py = []
    for y, det in enumerate(slices):
        Dc = K.ctx.from_dict({(0, 0): 1})
        for f, o in zip(facs, obs):
            if o:
                Dc = Dc * spec(f, vmin + y) ** o
        qq, r = divmod(Dc, det.d)
        assert r == 0, f'D* det not a polynomial on slice {y}'
        Py.append(det.n * qq)
    log(f'  D_c = D*; D* det is a polynomial on every slice [{time.time()-T0:.1f}s]')
    mons = sorted(set().union(*[set(p.to_dict().keys()) for p in Py]))
    V = flint.fmpq_mat(B + 2, B + 2, [flint.fmpq(y) ** l for y in range(B + 2) for l in range(B + 2)]).inv()
    Y = flint.fmpq_mat(B + 2, len(mons))
    for y, p in enumerate(Py):
        d = p.to_dict()
        for jx, mo in enumerate(mons):
            c = d.get(mo)
            if c is not None:
                Y[y, jx] = c
    Cf = V * Y
    above = sum(1 for jx in range(len(mons)) if Cf[B + 1, jx] != 0)
    assert above == 0, f'{above} nonzero coefficients of y^{B + 1}: the degree bound is violated'
    # N~ in the shifted variables: y (for v), x (for u = umin + x), s
    ctx3 = flint.fmpq_mpoly_ctx.get(('y', u, 's'), 'lex')
    terms = {}
    for jx, mo in enumerate(mons):
        for l in range(L + 1):
            c = Cf[l, jx]
            if c != 0:
                terms[(l,) + tuple(mo)] = c
    Nt = ctx3.from_dict(terms)
    g3 = ctx3.gens()
    Ns = Nt.compose(g3[0], g3[1] + C.smin[u], g3[2], ctx=ctx3)
    co = Ns.to_dict()
    nneg = sum(1 for c in co.values() if c < 0)
    base = any(c > 0 and e[0] == 0 and e[1] == 0 for e, c in co.items())
    log(f'  N~: {len(co)} terms after the shifts, degrees {Ns.degrees()}, {nneg} negative coefficients, a positive one '
        f'free of the sizes: {base}; no term in y^{B + 1} through the extra node [{time.time()-T0:.1f}s]')
    fpos = all(P.one_signed(f, C) == 1 for f, o in zip(facs, obs) if o)
    mp = C.minpoint()
    vals = [flint.fmpq(mp[w]) if w != 's' else flint.fmpq(1) for w in names2]
    A = [[(x.n(*vals) / x.d(*vals)) if x != 0 else flint.fmpq(0) for x in row] for row in Qg]
    pd = True
    for c in range(n):
        p = A[c][c]
        if p <= 0:
            pd = False
            break
        for r in range(c + 1, n):
            f = A[r][c] / p
            for cc in range(c + 1, n):
                A[r][cc] -= f * A[c][cc]
    log(f'  factors of D_c positive after the shifts: {fpos}; Qg positive definite at the minimal member, sigma = 1: '
        f'{pd} [{time.time()-T0:.1f}s]')
    ok = nneg == 0 and base and fpos and pd
    return ok, {'B': B, 'slice_var': v, 'L': L, 'Nterms': len(co), 'estar': estar, 'dexp': obs,
                'time': time.time() - T0, 'names': names}


if __name__ == '__main__':
    out = sys.argv[1]
    logf = lambda m: print(m, flush=True)
    T0 = time.time()
    results, infos = {}, {}
    allok = True
    for s in sys.argv[2:]:
        case, mode, fam = P.parse_fam(s)
        assert mode == 'sym' and sum(x is None for x in fam) == 2 and not (fam[2] is None and fam[3] is None)
        lab = P.fam_label(case, mode, fam)
        print(lab, flush=True)
        key = (case, mode) + tuple(('>=', T.SMIN[nm]) if x is None else x for x, nm in zip(fam, T.SIZES))
        try:
            ok, info = prove1(case, fam, logf)
            res = 'PD on complement (determinant route, one size sliced)' if ok else 'FAILED (a sign or the base point)'
        except AssertionError as ex:
            ok, info, res = False, {'time': 0}, f'FAILED ({ex})'
        allok &= ok
        results[key] = res
        infos[key] = info
        print(f'{lab} {res} ({info["time"]:.1f}s)', flush=True)
    with open(out, 'wb') as fh:
        pickle.dump({'results': results, 'info': infos, 'smin': dict(T.SMIN)}, fh)
    print(f'{sum(r.startswith("PD") for r in results.values())} families proved of {len(results)}; '
          f'{time.time() - T0:.0f}s')
    print('PASS' if allok else 'FAIL')
