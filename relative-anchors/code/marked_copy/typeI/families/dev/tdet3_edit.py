"""The determinant route for sym families with exactly three symbolic sizes that are not i and tau' together with a
third (so no structure lemma): two sizes sliced, the third and sigma free, with the irrep split (tsbuild.py).
NOTES.md, Section 8 (Proposition D3).

Q2g is the gauge-fixed first part of the split form and Q3 the second part.  Steps:
 prep   1. Q3 by greedy symmetric elimination with pivots one-signed after the shifts;
        2. D* = prod f^{e_f} over the irreducible factors f of the entry denominators of Q2g; T := D* det Q2g is a
           polynomial with deg_v T <= B_v = deg_v D* + max_p sum (deg_v num - deg_v den) for every size v;
        3. the sliced pair (v1, v2): the two sizes with the smallest (B_v1 + 1)(B_v2 + 1); u the third (free);
        4. D_c = prod f^{d_f}: d_f = e_f for the factors free of u and sigma, otherwise the largest multiplicity seen on
           a few slices, capped at e_f; R := D*/D_c is a polynomial; L_vk := B_vk - deg_vk R;
        5. Q2g positive definite at the minimal member and sigma = 1 (exact elimination).
 slices for y = 0..B_v1, z = 0..B_v2: det Q2g at v1 = min + y, v2 = min + z exactly over Q(u, sigma); D_c det must be
        a polynomial P_{y,z} there; saved to OUT/s_y_z.pkl (columns y claimed by the workers, saved slices skipped);
        a complete column is interpolated in z (stage A), and its coefficients of z^l with l > L_v2 must vanish.
 final  stage B: every coefficient interpolated in y through y = 0..B_v1; the coefficients of y^k with k > L_v1 must
        vanish.  So N~(y, z, u, sigma) has deg_y <= L_v1, deg_z <= L_v2 and N~ = P_{y,z} on the grid; T - R N~ has
        degree <= B_v1 in v1 and <= B_v2 in v2 and vanishes on the grid, so T = R N~ and det Q2g = N~/D_c.  Signs: after
        u = umin + x every coefficient of N~ is >= 0 and one free of y, z, x is > 0; every factor of D_c is one-signed
        positive after the shifts.  Then det Q2g > 0 on the connected real region, Q2g is continuous there, and with
        step 5 Q2g > 0 everywhere; with step 1 and tsbuild, Q >= 0 with kernel exactly span(kappa) at every member and
        every sigma > 0.
python tdet3.py prep FAM OUT | slices OUT W 0 | final OUT RESULT.pkl"""
import sys, os, time, pickle
import flint
from flintfield import Field
import tfamlib as T
import tprove as P
import tdet
import tsbuild

ORDER2 = {'MT': 'U0,U1,W0,W1', 'MU': 'U00,U01,U10,U11,W00,W01,W10,W11'}


def setup(case, fam, log, do_sanity):
    os.environ.setdefault('ORDER2', ORDER2[case])
    R = tsbuild.build_split(case, fam, log, do_sanity=do_sanity)
    C = R['C']
    n = len(R['keep2'])
    E = {(a, b): (x.n, x.d) for a, row in enumerate(R['Q2g']) for b, x in enumerate(row) if x != 0}
    facs, estar = tdet.factor_data(E, n)
    return R, C, E, C.syms, n, facs, estar


def spec_factor(f, names2, K, free, fix):
    gens = dict(zip(free, K.ctx.gens()))
    return f.compose(*[gens[w] if w in gens else K.ctx.from_dict({(0,) * len(free): fix[w]}) for w in names2],
                     ctx=K.ctx)


def famdir_meta(out):
    with open(os.path.join(out, 'meta.pkl'), 'rb') as fh:
        return pickle.load(fh)


def prep(spec, out, log=print):
    T0 = time.time()
    os.makedirs(out, exist_ok=True)
    case, mode, fam = P.parse_fam(spec)
    assert mode == 'sym' and sum(x is None for x in fam) == 3 and not (fam[2] is None and fam[3] is None)
    R, C, E, names2, n, facs, estar = setup(case, fam, log, True)
    n3 = len(R['three'])
    q3 = 'none'
    if n3:
        M = {(a, b): R['Q3'][a][b] for a in range(n3) for b in range(n3)}
        res, info = P.eliminate(M, list(range(n3)), C, time.time(), 1800, False, R['three'])
        assert res == 'OK', ('Q3 elimination', res)
        q3 = f'{n3} pivots one-signed: {info["pivots"]}'
    log(f'  Q3: {q3} [{time.time()-T0:.1f}s]')
    szs = [v for v in names2 if v != 's']
    B = {v: tdet.degree_bound(E, n, facs, estar, names2.index(v)) for v in szs}
    pairs = [(a, b) for a in szs for b in szs if a != b]
    v1, v2 = min(pairs, key=lambda p: ((B[p[0]] + 1) * (B[p[1]] + 1), -B[p[0]]))
    u = [v for v in szs if v not in (v1, v2)][0]
    log(f'  Q2g {n} x {n}; {len(facs)} irreducible factors, exponents of D* {estar}; bounds {B}; sliced {v1}, {v2}; '
        f'free {u} [{time.time()-T0:.1f}s]')
    free = [u, 's']
    K = Field(free)
    const = [f.degrees()[names2.index(u)] == 0 and f.degrees()[names2.index('s')] == 0 for f in facs]
    B1, B2 = B[v1], B[v2]
    probe = sorted({(B1 // 3, B2 // 2), (B1, B2), (B1 // 2, B2 - 1), (2 * B1 // 3, B2 // 3)})
    obs = [0] * len(facs)
    for (y, z) in probe:
        fix = {v1: C.smin[v1] + y, v2: C.smin[v2] + z}
        K2, det = tdet.slice_det(E, n, None, names2, fix, free)
        for k_, f in enumerate(facs):
            if not const[k_]:
                obs[k_] = max(obs[k_], tdet.multiplicity(spec_factor(f, names2, K, free, fix), det.d))
    dexp = [e if c else min(o, e) for e, c, o in zip(estar, const, obs)]
    degR = {v: sum((e - d) * f.degrees()[names2.index(v)] for f, e, d in zip(facs, estar, dexp)) for v in (v1, v2)}
    L1, L2 = B1 - degR[v1], B2 - degR[v2]
    log(f'  observed multiplicities on the probes {probe}: {obs}; D_c exponents {dexp}; deg R {degR}; '
        f'L = ({L1}, {L2}) [{time.time()-T0:.1f}s]')
    mp = C.minpoint()
    vals = [flint.fmpq(mp[w]) if w != 's' else flint.fmpq(1) for w in C.syms]
    A = [[(x.n(*vals) / x.d(*vals)) if x != 0 else flint.fmpq(0) for x in row] for row in R['Q2g']]
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
    assert pd, 'Q2g not positive definite at the minimal member, sigma = 1'
    fpos = all(P.one_signed(f, C) == 1 for f, d in zip(facs, dexp) if d)
    assert fpos, 'a factor of D_c is not one-signed positive after the shifts'
    meta = {'spec': spec, 'case': case, 'fam': fam, 'v1': v1, 'v2': v2, 'u': u, 'B1': B1, 'B2': B2, 'L1': L1,
            'L2': L2, 'degR': degR, 'facs': [str(f) for f in facs], 'estar': estar, 'dexp': dexp,
            'names': R['names'], 'keep2': R['keep2'], 'three': R['three'], 'gauge': R['gauge'], 'q3': q3,
            'npts': R['npts'], 'base_pd': pd, 'fpos': fpos, 'bounds': B, 'prep_time': time.time() - T0}
    with open(os.path.join(out, 'meta.pkl'), 'wb') as fh:
        pickle.dump(meta, fh)
    log(f'  Q2g positive definite at the minimal member, sigma = 1; factors of D_c positive after the shifts; '
        f'grid {B1 + 1} x {B2 + 1} = {(B1 + 1) * (B2 + 1)} slices [{time.time()-T0:.1f}s]')
    return meta


def inv_vandermonde(D):
    V = flint.fmpq_mat(D + 1, D + 1, [flint.fmpq(z) ** l for z in range(D + 1) for l in range(D + 1)])
    return V.inv()


def stage_a(out, y, B2, L2, log=print):
    t0 = time.time()
    cols, mons = [], set()
    for z in range(B2 + 1):
        with open(os.path.join(out, f's_{y}_{z}.pkl'), 'rb') as fh:
            d = pickle.load(fh)
        cols.append(d)
        mons.update(d.keys())
    mons = sorted(mons)
    Y = flint.fmpq_mat(B2 + 1, len(mons))
    for z, d in enumerate(cols):
        for j, mo in enumerate(mons):
            if mo in d:
                p, q = d[mo]
                Y[z, j] = flint.fmpq(p, q)
    Cz = inv_vandermonde(B2) * Y
    data, bad = {}, 0
    for j, mo in enumerate(mons):
        for l in range(B2 + 1):
            c = Cz[l, j]
            if c != 0:
                if l > L2:
                    bad += 1
                else:
                    data[(l,) + tuple(mo)] = (int(c.p), int(c.q))
    fn = os.path.join(out, f'col_{y}.pkl')
    with open(fn + '.tmp', 'wb') as fh:
        pickle.dump({'bad': bad, 'data': data}, fh, protocol=4)
    os.replace(fn + '.tmp', fn)
    log(f'  column y={y}: {len(mons)} monomials, {len(data)} coefficients, {bad} nonzero above z^{L2} '
        f'[{time.time()-t0:.1f}s]')
    return bad == 0


def slices(out, w, log=print):
    T0 = time.time()
    meta = famdir_meta(out)
    R, C, E, names2, n, facs, estar = setup(meta['case'], meta['fam'], lambda m: None, False)
    assert [str(f) for f in facs] == meta['facs'] and estar == meta['estar'] and R['keep2'] == meta['keep2']
    v1, v2, u, B1, B2, L2 = meta['v1'], meta['v2'], meta['u'], meta['B1'], meta['B2'], meta['L2']
    free = [u, 's']
    K = Field(free)
    ok, nd = True, 0
    for y in range(B1 + 1):
        if os.path.exists(os.path.join(out, 'STOP')) or os.path.exists('STOP_ALL'):
            log(f'{meta["spec"]}: stop file found, no further columns claimed')
            break
        try:
            fd = os.open(os.path.join(out, f'claim_{y}'), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.write(fd, str(w).encode())
            os.close(fd)
        except FileExistsError:
            continue
        for z in range(B2 + 1):
            fn = os.path.join(out, f's_{y}_{z}.pkl')
            if os.path.exists(fn):
                continue
            t1 = time.time()
            fix = {v1: C.smin[v1] + y, v2: C.smin[v2] + z}
            K2, det = tdet.slice_det(E, n, None, names2, fix, free)
            Dc = K.ctx.from_dict({(0, 0): 1})
            for f, d in zip(facs, meta['dexp']):
                if d:
                    Dc = Dc * spec_factor(f, names2, K, free, fix) ** d
            q, r = divmod(Dc, det.d)
            if r != 0 or det == 0:
                log(f'{meta["spec"]} y={y} z={z}: D_c det NOT A POLYNOMIAL (or det = 0)')
                ok = False
                continue
            Pl = det.n * q
            data = {tuple(e): (int(c.p), int(c.q)) for e, c in Pl.to_dict().items()}
            with open(fn + '.tmp%d' % w, 'wb') as fh:
                pickle.dump(data, fh, protocol=4)
            try:
                os.replace(fn + '.tmp%d' % w, fn)
            except OSError:
                if not os.path.exists(fn):
                    raise
                os.remove(fn + '.tmp%d' % w)
            nd += 1
            if z == B2 or z == 0:
                log(f'{meta["spec"]} y={y} z={z}: {len(data)} terms [{time.time()-t1:.1f}s, total {time.time()-T0:.0f}s]')
        if all(os.path.exists(os.path.join(out, f's_{y}_{z}.pkl')) for z in range(B2 + 1)) and \
                not os.path.exists(os.path.join(out, f'col_{y}.pkl')):
            ok &= stage_a(out, y, B2, L2, log)
    log(f'{meta["spec"]} worker {w}: {nd} slices, {time.time()-T0:.0f}s')
    return ok


def final(out, result, log=print, chunk=4000):
    T0 = time.time()
    meta = famdir_meta(out)
    done = os.path.join(out, 'final_ok')
    if os.path.exists(done) and os.path.exists(result) and os.environ.get('REDO_FINAL') != '1':
        # a final step that already passed for this family (marker written below); REDO_FINAL=1 recomputes it
        log(open(done).read().strip())
        log('(final step already done; set REDO_FINAL=1 to recompute)')
        return True
    spec, B1, B2, L1, L2, u = meta['spec'], meta['B1'], meta['B2'], meta['L1'], meta['L2'], meta['u']
    case, fam = meta['case'], meta['fam']
    R, C, E, names2, n, facs, estar = setup(case, fam, log, True)
    assert [str(f) for f in facs] == meta['facs'] and estar == meta['estar'] and R['keep2'] == meta['keep2']
    missing = [(y, z) for y in range(B1 + 1) for z in range(B2 + 1)
               if not os.path.exists(os.path.join(out, f's_{y}_{z}.pkl'))]
    assert not missing, f'{len(missing)} slices missing, e.g. {missing[:5]}'
    abad = 0
    for y in range(B1 + 1):
        if not os.path.exists(os.path.join(out, f'col_{y}.pkl')):
            stage_a(out, y, B2, L2, log)
        with open(os.path.join(out, f'col_{y}.pkl'), 'rb') as fh:
            abad += pickle.load(fh)['bad']
    log(f'  all {(B1 + 1) * (B2 + 1)} slices and {B1 + 1} columns present; nonzero coefficients above z^{L2}: {abad} '
        f'[{time.time()-T0:.0f}s]')
    Vy = inv_vandermonde(B1)
    keys = set()
    for y in range(B1 + 1):
        with open(os.path.join(out, f'col_{y}.pkl'), 'rb') as fh:
            keys.update(pickle.load(fh)['data'].keys())
    keys = sorted(keys)
    ctx = flint.fmpq_mpoly_ctx.get(('y', 'z', u, 's'), 'lex')
    gy, gz, gu, gs = ctx.gens()
    nbad = nneg = npos = 0
    base = False
    degs = [0, 0, 0, 0]
    umin = C.smin[u]
    # chunks: all keys with one z-power l (so every group (k, l) below holds all powers of u and s of N~ for that
    # (y-power, z-power), and the shift u = umin + x is applied to complete polynomials)
    byl = {}
    for key in keys:
        byl.setdefault(key[0], []).append(key)
    for c0, lz in enumerate(sorted(byl)):
        part = byl[lz]
        Y = flint.fmpq_mat(B1 + 1, len(part))
        for y in range(B1 + 1):
            with open(os.path.join(out, f'col_{y}.pkl'), 'rb') as fh:
                data = pickle.load(fh)['data']
            for j, mo in enumerate(part):
                pq = data.get(mo)
                if pq is not None:
                    Y[y, j] = flint.fmpq(*pq)
            del data
        Cy = Vy * Y
        # the terms of this chunk, grouped by (y-power, z-power): each a polynomial in (u, s), shifted u = umin + x
        groups = {}
        for j, (l, a, b) in enumerate(part):
            for k in range(B1 + 1):
                c = Cy[k, j]
                if c == 0:
                    continue
                if k > L1:
                    nbad += 1
                    continue
                groups.setdefault((k, l), {})[(a, b)] = c
        del Y, Cy
        ctx2 = flint.fmpq_mpoly_ctx.get((u, 's'), 'lex')
        xx, ss = ctx2.gens()
        for (k, l), g in groups.items():
            Pa = ctx2.from_dict(g)
            Pb = Pa.compose(xx + umin, ss, ctx=ctx2)
            for (a, b), c in Pb.to_dict().items():
                if c < 0:
                    nneg += 1
                elif c > 0:
                    npos += 1
                    degs = [max(degs[0], k), max(degs[1], l), max(degs[2], a), max(degs[3], b)]
                    if k == 0 and l == 0 and a == 0:
                        base = True
        log(f'    z^{lz}: {len(part)} monomials; so far nonzero above y^{L1} {nbad}, positive {npos}, negative {nneg} '
            f'[{time.time()-T0:.0f}s]')
    log(f'  N~(y, z, x, s) after {u} = {umin} + x: {npos} positive and {nneg} negative coefficients, degrees <= '
        f'{tuple(degs)}, a positive one free of y, z, x: {base}; nonzero above y^{L1}: {nbad} [{time.time()-T0:.0f}s]')
    ok = abad == 0 and nbad == 0 and nneg == 0 and base and meta['base_pd'] and meta['fpos']
    key = (case, 'sym') + tuple(('>=', T.SMIN[nm]) if x is None else x for x, nm in zip(fam, T.SIZES))
    res = 'PD on complement (determinant route with the split, two sizes sliced, no structure)' if ok else 'FAILED'
    info = dict(meta)
    info.update({'Nterms_shifted': npos + nneg, 'Ndeg': degs, 'time_final': time.time() - T0})
    with open(result, 'wb') as fh:
        pickle.dump({'results': {key: res}, 'info': {key: info}, 'smin': dict(T.SMIN)}, fh)
    log(f'{P.fam_label(case, "sym", fam)} {res}')
    if ok:
        with open(done, 'w') as fh:
            fh.write(f'{P.fam_label(case, "sym", fam)} {res}: {npos} positive and {nneg} negative coefficients '
                     f'after the shift, degrees {tuple(degs)}, nonzero above y^{L1}: {nbad}, above z^{L2}: {abad}'
                     + chr(10))
    return ok


if __name__ == '__main__':
    logf = lambda m: print(m, flush=True)
    mode = sys.argv[1]
    if mode == 'prep':
        print(sys.argv[2], flush=True)
        prep(sys.argv[2], sys.argv[3], logf)
        print('PASS')
    elif mode == 'slices':
        ok = slices(sys.argv[2], int(sys.argv[3]), logf)
        print('PASS' if ok else 'FAIL')
    elif mode == 'final':
        ok = final(sys.argv[2], sys.argv[3], logf)
        print('PASS' if ok else 'FAIL')
