"""The determinant route for sym families with i and tau' symbolic and exactly one of q0, t symbolic (call it v), with
the irrep split (tsbuild.py).  NOTES.md, Section 8 (Proposition D2).

Q2g is the gauge-fixed first part of the split form and Q3 the second part (tsbuild.build_split: preliminaries,
block-diagonality, kernel vectors).  Steps:
 prep   1. Q3 by greedy symmetric elimination with pivots one-signed after the shifts (tprove.eliminate);
        2. structure (Lemma N): with tau' = nu - 2i every entry of Q2g is p/d, d free of i, deg_i p <= 2 (tdet.to_nu);
        3. D* = prod f^{e_f} over the irreducible factors f of the entry denominators (e_f: the largest sum of
           multiplicities along a permutation); T := D* det Q2g is a polynomial with deg_i T <= B_i (largest sum of
           deg_i p along a permutation) and deg_v T <= B_v = deg_v D* + max_p sum (deg_v num - deg_v den);
        4. D_c = prod f^{d_f}: d_f = e_f for the factors free of nu and sigma (constant on every slice), otherwise the
           largest multiplicity seen on a few slices, capped at e_f, so R := D*/D_c is a polynomial, free of i;
           L_v := B_v - deg_v R;
        5. Q2g positive definite at the minimal member and sigma = 1 (exact elimination).
 slices for y = 0..B_v, z = 0..B_i: det Q2g at v = vmin + y, i = 2 + z exactly over Q(nu, sigma); D_c det must be a
        polynomial P_{y,z} there (exact division); saved to OUT/s_y_z.pkl.  Worker W of NW takes the columns y with
        y % NW == W (NW = 0: each worker claims the next free column); a slice already saved is skipped; a complete
        column is interpolated in z (stage A).
 final  stage A for every column (interpolant of degree <= B_i in z through z = 0..B_i, inverse Vandermonde over Q);
        stage B: every coefficient interpolated in y through y = 0..B_v; the coefficients of y^k with k > L_v must
        vanish.  The result N~(y, z, nu, sigma) has deg_y <= L_v, deg_z <= B_i and N~ = P_{y,z} on the grid.  Then
        T - R N~ (deg_v <= B_v since deg_v R + L_v = B_v, deg_i <= B_i) vanishes on the grid, so T = R N~ and
        det Q2g = N~ / D_c identically.  Signs: with nu = tau' + 2i, i = 2 + z, tau' = 2 + w, every coefficient of
        N~(y, z, w, sigma) is >= 0 and one free of y, z, w is > 0; every factor of D_c is one-signed positive after the
        shifts.  So det Q2g > 0 on the connected real region, Q2g is continuous there, and with step 5 Q2g > 0
        on the whole region; with step 1 and tsbuild, Q >= 0 with kernel exactly span(kappa) at every member and
        every sigma > 0.
python tdet2.py prep FAM OUT | slices OUT W NW | final OUT RESULT.pkl"""
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
    E, ctx2, names2 = tdet.to_nu(C, R['Q2g'])
    n = len(R['keep2'])
    facs, estar = tdet.factor_data(E, n)
    return R, C, E, ctx2, names2, n, facs, estar


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
    assert mode == 'sym' and fam[2] is None and fam[3] is None and sum(x is None for x in fam[:2]) == 1
    R, C, E, ctx2, names2, n, facs, estar = setup(case, fam, log, True)
    v = 'q0' if fam[0] is None else 't'
    # 1. Q3 by elimination
    n3 = len(R['three'])
    q3 = 'none'
    if n3:
        M = {(a, b): R['Q3'][a][b] for a in range(n3) for b in range(n3)}
        res, info = P.eliminate(M, list(range(n3)), C, time.time(), 1800, False, R['three'])
        assert res == 'OK', ('Q3 elimination', res)
        q3 = f'{n3} pivots one-signed: {info["pivots"]}'
    log(f'  Q3: {q3} [{time.time()-T0:.1f}s]')
    # 3. bounds
    ii = names2.index('i')
    Wi = [[(E[(a, b)][0].degrees()[ii] if (a, b) in E else None) for b in range(n)] for a in range(n)]
    Bi = tdet.maxperm(Wi)
    iv = names2.index(v)
    Bv = tdet.degree_bound(E, n, facs, estar, iv)
    log(f'  Q2g {n} x {n}; {len(facs)} irreducible factors, exponents of D* {estar}; deg_i T <= B_i = {Bi}, '
        f'deg_{v} T <= B_v = {Bv} [{time.time()-T0:.1f}s]')
    # 4. observed exponents on a few slices
    free = ['nu', 's']
    K = Field(free)
    vmin = C.smin[v]
    obs = [0] * len(facs)
    const = []
    for f in facs:
        dg = f.degrees()
        const.append(dg[names2.index('nu')] == 0 and dg[names2.index('s')] == 0)
    # probes at members away from the smallest sizes (where specialised factors may coincide); d_f is capped at e_f,
    # and every slice checks that D_c det is a polynomial
    probe = sorted({(Bv // 3, Bi // 2), (Bv, Bi), (Bv // 2, Bi - 1), (2 * Bv // 3, Bi // 3)})
    for (y, z) in probe:
        fix = {v: vmin + y, 'i': 2 + z}
        K2, det = tdet.slice_det(E, n, ctx2, names2, fix, free)
        for k_, f in enumerate(facs):
            if const[k_]:
                continue
            fs = spec_factor(f, names2, K, free, fix)
            obs[k_] = max(obs[k_], tdet.multiplicity(fs, det.d))
    dexp = [e if c else min(o, e) for e, c, o in zip(estar, const, obs)]
    log(f'  observed multiplicities on the probes {probe}: {obs}')
    degR = sum((e - d) * f.degrees()[iv] for f, e, d in zip(facs, estar, dexp))
    Lv = Bv - degR
    log(f'  D_c exponents {dexp} (constant on slices: {[k for k, c in enumerate(const) if c]}); deg_{v} R = {degR}, '
        f'L_v = {Lv} [{time.time()-T0:.1f}s]')
    # 5. base point
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
    meta = {'spec': spec, 'case': case, 'fam': fam, 'v': v, 'vmin': vmin, 'Bv': Bv, 'Bi': Bi, 'Lv': Lv,
            'degR': degR, 'facs': [str(f) for f in facs], 'estar': estar, 'dexp': dexp, 'names': R['names'],
            'keep2': R['keep2'], 'three': R['three'], 'gauge': R['gauge'], 'q3': q3, 'npts': R['npts'],
            'base_pd': pd, 'fpos': fpos, 'prep_time': time.time() - T0}
    with open(os.path.join(out, 'meta.pkl'), 'wb') as fh:
        pickle.dump(meta, fh)
    log(f'  Q2g positive definite at the minimal member, sigma = 1; factors of D_c positive after the shifts; '
        f'grid {Bv + 1} x {Bi + 1} = {(Bv + 1) * (Bi + 1)} slices [{time.time()-T0:.1f}s]')
    return meta


def slices(out, w, nw, log=print):
    T0 = time.time()
    meta = famdir_meta(out)
    R, C, E, ctx2, names2, n, facs, estar = setup(meta['case'], meta['fam'], lambda m: None, False)
    assert [str(f) for f in facs] == meta['facs'] and estar == meta['estar'] and R['keep2'] == meta['keep2']
    v, vmin, Bv, Bi = meta['v'], meta['vmin'], meta['Bv'], meta['Bi']
    free = ['nu', 's']
    K = Field(free)
    ok = True
    nd = 0
    for y in range(Bv + 1):
        if nw > 0 and y % nw != w:
            continue
        if os.path.exists(os.path.join(out, 'STOP')) or os.path.exists('STOP_ALL'):
            log(f'{meta["spec"]}: stop file found, no further columns claimed')
            break
        if nw == 0:
            # dynamic assignment: claim the column y by creating OUT/claim_y exclusively (run scripts remove the claims
            # of incomplete columns before starting workers, so a stopped run resumes cleanly)
            try:
                fd = os.open(os.path.join(out, f'claim_{y}'), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                os.write(fd, str(w).encode())
                os.close(fd)
            except FileExistsError:
                continue
        for z in range(Bi + 1):
            fn = os.path.join(out, f's_{y}_{z}.pkl')
            if os.path.exists(fn):
                continue
            t1 = time.time()
            fix = {v: vmin + y, 'i': 2 + z}
            K2, det = tdet.slice_det(E, n, ctx2, names2, fix, free)
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
            if z == Bi or z == 0:
                log(f'{meta["spec"]} y={y} z={z}: {len(data)} terms [{time.time()-t1:.1f}s, total {time.time()-T0:.0f}s]')
        if all(os.path.exists(os.path.join(out, f's_{y}_{z}.pkl')) for z in range(Bi + 1)) and \
                not os.path.exists(os.path.join(out, f'col_{y}.pkl')):
            ok &= stage_a(out, y, Bi, log)
    log(f'{meta["spec"]} worker {w}/{nw}: {nd} slices, {time.time()-T0:.0f}s')
    return ok


def inv_vandermonde(D):
    V = flint.fmpq_mat(D + 1, D + 1, [flint.fmpq(z) ** l for z in range(D + 1) for l in range(D + 1)])
    return V.inv()


def stage_a(out, y, Bi, log=print):
    t0 = time.time()
    cols, mons = [], set()
    for z in range(Bi + 1):
        with open(os.path.join(out, f's_{y}_{z}.pkl'), 'rb') as fh:
            d = pickle.load(fh)
        cols.append(d)
        mons.update(d.keys())
    mons = sorted(mons)
    Y = flint.fmpq_mat(Bi + 1, len(mons))
    for z, d in enumerate(cols):
        for j, mo in enumerate(mons):
            if mo in d:
                p, q = d[mo]
                Y[z, j] = flint.fmpq(p, q)
    Cz = inv_vandermonde(Bi) * Y
    data = {}
    for j, mo in enumerate(mons):
        for l in range(Bi + 1):
            c = Cz[l, j]
            if c != 0:
                data[(l,) + tuple(mo)] = (int(c.p), int(c.q))
    fn = os.path.join(out, f'col_{y}.pkl')
    with open(fn + '.tmp', 'wb') as fh:
        pickle.dump(data, fh, protocol=4)
    os.replace(fn + '.tmp', fn)
    log(f'  column y={y}: {len(mons)} monomials in (nu, s), {len(data)} coefficients after interpolation in z '
        f'[{time.time()-t0:.1f}s]')
    return True


def final(out, result, log=print, chunk=8000):
    T0 = time.time()
    meta = famdir_meta(out)
    done = os.path.join(out, 'final_ok')
    if os.path.exists(done) and os.path.exists(result) and os.environ.get('REDO_FINAL') != '1':
        # a final step that already passed for this family (marker written below); REDO_FINAL=1 recomputes it
        log(open(done).read().strip())
        log('(final step already done; set REDO_FINAL=1 to recompute)')
        return True
    spec, Bv, Bi, Lv = meta['spec'], meta['Bv'], meta['Bi'], meta['Lv']
    case, fam = meta['case'], meta['fam']
    # recompute the preliminaries from scratch (sanity included) and compare with the meta data
    R, C, E, ctx2, names2, n, facs, estar = setup(case, fam, log, True)
    assert [str(f) for f in facs] == meta['facs'] and estar == meta['estar'] and R['keep2'] == meta['keep2']
    missing = [(y, z) for y in range(Bv + 1) for z in range(Bi + 1)
               if not os.path.exists(os.path.join(out, f's_{y}_{z}.pkl'))]
    assert not missing, f'{len(missing)} slices missing, e.g. {missing[:5]}'
    for y in range(Bv + 1):
        if not os.path.exists(os.path.join(out, f'col_{y}.pkl')):
            stage_a(out, y, Bi, log)
    log(f'  all {(Bv + 1) * (Bi + 1)} slices and {Bv + 1} columns present [{time.time()-T0:.0f}s]')
    Vy = inv_vandermonde(Bv)
    keys = set()
    for y in range(Bv + 1):
        with open(os.path.join(out, f'col_{y}.pkl'), 'rb') as fh:
            keys.update(pickle.load(fh).keys())
    keys = sorted(keys)
    log(f'  {len(keys)} monomials z^l nu^a s^b over the columns [{time.time()-T0:.0f}s]')
    # stage B, in chunks of monomials made of whole groups with one power b of s: the substitution nu = 6 + 2z + w
    # mixes the powers of z and nu, so it is applied, for each (y-power k, s-power b), to the complete polynomial in
    # (z, nu), and every such polynomial lies in one chunk
    byb = {}
    for key in keys:
        byb.setdefault(key[2], []).append(key)
    chunks, cur = [], []
    for b in sorted(byb):
        if cur and len(cur) + len(byb[b]) > chunk:
            chunks.append(cur)
            cur = []
        cur = cur + byb[b]
    if cur:
        chunks.append(cur)
    ctxA = flint.fmpq_mpoly_ctx.get(('z', 'nu'), 'lex')
    ctxB = flint.fmpq_mpoly_ctx.get(('z', 'w'), 'lex')
    zz, ww = ctxB.gens()
    nbad = nterms = nneg = npos = 0
    base = False
    degs = [0, 0, 0, 0]
    for ci, part in enumerate(chunks):
        Y = flint.fmpq_mat(Bv + 1, len(part))
        for y in range(Bv + 1):
            with open(os.path.join(out, f'col_{y}.pkl'), 'rb') as fh:
                data = pickle.load(fh)
            for j, mo in enumerate(part):
                pq = data.get(mo)
                if pq is not None:
                    Y[y, j] = flint.fmpq(*pq)
            del data
        Cy = Vy * Y
        del Y
        groups = {}
        for j, (l, a, b) in enumerate(part):
            for k in range(Bv + 1):
                c = Cy[k, j]
                if c == 0:
                    continue
                if k > Lv:
                    nbad += 1
                    continue
                groups.setdefault((k, b), {})[(l, a)] = c
        del Cy
        for (k, b), g in groups.items():
            nterms += len(g)
            Pa = ctxA.from_dict(g)
            Pb = Pa.compose(zz, 6 + 2 * zz + ww, ctx=ctxB)
            for (l, ww_), c in Pb.to_dict().items():
                if c < 0:
                    nneg += 1
                elif c > 0:
                    npos += 1
                    degs = [max(degs[0], k), max(degs[1], l), max(degs[2], ww_), max(degs[3], b)]
                    if k == 0 and l == 0 and ww_ == 0:
                        base = True
        del groups
        log(f'    chunk {ci + 1}/{len(chunks)} ({len(part)} monomials) interpolated in y; so far nonzero above '
            f'y^{Lv}: {nbad}, after the substitution {npos} positive, {nneg} negative [{time.time()-T0:.0f}s]')
    log(f'  N~(y, z, nu, s): {nterms} terms; coefficients of y^k with k > L_v = {Lv}: {nbad} nonzero '
        f'[{time.time()-T0:.0f}s]')
    log(f'  N~ with nu = tau1 + 2i, i = 2 + z, tau1 = 2 + w: {npos} positive and {nneg} negative coefficients, '
        f'degrees (y, z, w, s) <= {tuple(degs)}, a positive one free of y, z, w: {base} [{time.time()-T0:.0f}s]')
    ok = nbad == 0 and nneg == 0 and base and meta['base_pd'] and meta['fpos']
    key = (case, 'sym') + tuple(('>=', T.SMIN[nm]) if x is None else x for x, nm in zip(fam, T.SIZES))
    res = 'PD on complement (determinant route with the split, two sizes sliced)' if ok else 'FAILED'
    info = dict(meta)
    info.update({'Nterms_shifted': npos + nneg, 'Ndeg': degs, 'time_final': time.time() - T0})
    with open(result, 'wb') as fh:
        pickle.dump({'results': {key: res}, 'info': {key: info}, 'smin': dict(T.SMIN)}, fh)
    log(f'{P.fam_label(case, "sym", fam)} {res}')
    if ok:
        with open(done, 'w') as fh:
            fh.write(f'{P.fam_label(case, "sym", fam)} {res}: {npos} positive and {nneg} negative coefficients '
                     f'after the substitution, degrees {tuple(degs)}, nonzero above y^{Lv}: {nbad}' + chr(10))
    return ok


if __name__ == '__main__':
    logf = lambda m: print(m, flush=True)
    mode = sys.argv[1]
    if mode == 'prep':
        print(sys.argv[2], flush=True)
        prep(sys.argv[2], sys.argv[3], logf)
        print('PASS')
    elif mode == 'slices':
        ok = slices(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), logf)
        print('PASS' if ok else 'FAIL')
    elif mode == 'final':
        ok = final(sys.argv[2], sys.argv[3], logf)
        print('PASS' if ok else 'FAIL')
