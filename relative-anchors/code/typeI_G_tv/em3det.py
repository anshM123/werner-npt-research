"""Theorem TV for both level sizes symbolic: the third LDL^T pivot of (B_inv) by the determinant route
(NOTES.md, Section 10).  A is the 3 x 3 matrix of (B_inv) (em3lib.build_A; u_i = 2 + x_i, sigma = s_t + r, s1 = s),
for s_t = s (KA = s) and s_t = 1 (KA = 1).  Region: x1, x2 >= 0, r > 0, s > 0.

 bounds KA   every entry denominator of A is a product of irreducible factors with positive coefficients and a monomial
             free of x1, x2, r (so A is continuous on the region); the exponents e_f of D*, the bound B for deg_x1 T,
             T = D* det A; the exponents d_f of the candidate denominator Den_c (seen on the slices x1 = 0, 7, 31) with
             d_f <= e_f; deg_x1 R for R = D*/Den_c and the low degree L = B - deg_x1 R; every factor of Den_c is
             positive on the region.  Written to em3_out/setup_KA.pkl.
 slices KA   for c = 0..B: det A at x1 = c exactly over Q(x2, r, s); Den_c(x1 = c) * det must be a polynomial P_c
             (exact division); P_c is written to em3_out/KA/s_c.pkl.
 interp KA   for every monomial x2^a r^b s^e: the values of its coefficient in P_c, c = 0..B, are interpolated exactly in
             x1 (inverse Vandermonde over Q, nodes 0..B); every coefficient of x1^l with l > L must vanish.  The result
             N~(x1, x2, r, s) (degree <= L in x1) must have only coefficients >= 0 and a positive coefficient free of
             x1, x2, r.  Written to em3_out/N_KA.pkl.
 final       both KA: the setup recomputed from scratch and compared, the first two pivots of A certified again from
             the same matrix (prove_gersh.cert), every slice present, the interpolation passed.
Conclusion: T and N~ R are polynomials of degree <= B in x1 that agree at x1 = 0..B (as polynomials in x2, r, s), so
T = N~ R and det A = N~/Den_c > 0 on the region; with the first two pivots positive, the third pivot is positive.
python em3det.py bounds|slices|interp KA      python em3det.py final"""
import sys, os, time, pickle, gc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
import em3lib as L

Q = flint.fmpq
OUT = 'em3_out'
FREE = ('x1', 'x2', 'r')
OBS_SLICES = (0, 7, 31)


def ser(p):
    return {tuple(int(e) for e in m): (int(c.p), int(c.q)) for m, c in zip(p.monoms(), p.coeffs())}


def setup(ka, log=print):
    """the bookkeeping of the mode bounds, from scratch; returns (A, Ao-slicer, data)"""
    t0 = time.time()
    A, (names, Ao) = L.build_A(ka, orig=True, log=log)
    OS = L.OrigSlicer(ka, names, Ao)
    factors, mult, estar = L.factor_data(A)
    dens_ok = all(L.one_signed_pos(f, FREE) for f in factors)
    bounds = {v: L.degree_bound(A, factors, estar, iv) for iv, v in enumerate(L.VN)}
    log('  %d irreducible factors of the entry denominators, all with positive coefficients and a monomial free of '
        'x1, x2, r: %s; exponents of D*: %s [%.1fs]' % (len(factors), dens_ok, estar, time.time() - t0), flush=True)
    data = dict(ka=ka, facs=[str(f) for f in factors], fac_dicts=[ser(f) for f in factors], estar=estar,
                dens_ok=dens_ok, bounds={v: b for v, (b, d) in bounds.items()},
                degDstar={v: d for v, (b, d) in bounds.items()},
                mult={k: v for k, v in mult.items()})
    return A, OS, factors, data


def den_c(factors, dexp):
    D = L.ONE
    for f, e in zip(factors, dexp):
        D = D * f ** e
    return D


def mode_bounds(ka):
    t0 = time.time()
    A, OS, factors, data = setup(ka)
    obs = [0] * len(factors)
    for c in OBS_SLICES:
        t1 = time.time()
        dt = OS.det(c)
        row = []
        for k, f in enumerate(factors):
            fs = L.spec(f, c)
            if fs.is_constant():
                row.append('-')
                continue
            m = L.multiplicity(fs, dt.d)
            assert L.multiplicity(fs, dt.n) == 0
            obs[k] = max(obs[k], m)
            row.append(str(m))
        print('  slice x1 = %d: det numerator %d terms, denominator %d terms; exponents of the factors in the '
              'denominator: %s [%.1fs]' % (c, len(dt.n), len(dt.d), ' '.join(row), time.time() - t1), flush=True)
    estar = data['estar']
    divides = all(d <= e for d, e in zip(obs, estar))
    B = data['bounds']['x1']
    dDs = data['degDstar']['x1']
    dDc = sum(e * f.degrees()[0] for f, e in zip(factors, obs))
    degR = dDs - dDc
    low = B - degR
    dc_ok = all(L.one_signed_pos(f, FREE) for f, e in zip(factors, obs) if e > 0)
    print('  exponents of D*: %s; of Den_c: %s; Den_c divides D*: %s' % (estar, obs, divides))
    print('  degree bounds of T = D* det A (x1, x2, r, s): %s; deg D*: %s' % (data['bounds'], data['degDstar']))
    print('  x1: B = %d, deg D* = %d, deg Den_c = %d, deg R = %d, low degree L = B - deg R = %d; grid x1 = 0..%d '
          '(%d slices)' % (B, dDs, dDc, degR, low, B, B + 1))
    print('  every factor of Den_c has positive coefficients and a monomial free of x1, x2, r: %s' % dc_ok)
    data.update(dexp=obs, B=B, degR=degR, L=low, divides=divides, dc_ok=dc_ok)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, 'setup_%s.pkl' % ka), 'wb') as fh:
        pickle.dump(data, fh)
    ok = data['dens_ok'] and divides and dc_ok
    print('bounds s_t = %s [%.1fs]' % (ka, time.time() - t0))
    return ok


def load_setup(ka):
    with open(os.path.join(OUT, 'setup_%s.pkl' % ka), 'rb') as fh:
        return pickle.load(fh)


def mode_slices(ka, resume=False):
    t0 = time.time()
    st = load_setup(ka)
    A, OS, factors, data = setup(ka, log=lambda *a, **k: None)
    assert data['facs'] == st['facs'] and data['estar'] == st['estar']
    B = st['B']
    d = os.path.join(OUT, ka)
    os.makedirs(d, exist_ok=True)
    good = True
    npos_all = True
    reused = 0
    for c in range(B + 1):
        fn = os.path.join(d, 's_%d.pkl' % c)
        if resume and os.path.exists(fn):
            reused += 1
            continue
        t1 = time.time()
        dt = OS.det(c)
        Dc = L.ONE3
        for f, e in zip(factors, st['dexp']):
            if e:
                Dc = Dc * L.spec(f, c) ** e
        q, rem = divmod(Dc, dt.d)
        if rem != 0:
            print('x1 = %d: Den_c * det is NOT a polynomial' % c, flush=True)
            good = False
            continue
        P = dt.n * q
        cf = list(P.coeffs())
        allpos = all(x > 0 for x in cf)
        npos_all &= allpos
        with open(fn + '.tmp', 'wb') as fh:
            pickle.dump(ser(P), fh, protocol=4)
        os.replace(fn + '.tmp', fn)
        print('x1 = %d: Den_c det is a polynomial, %d terms, all coefficients positive: %s [%.1fs, total %.0fs]'
              % (c, len(cf), allpos, time.time() - t1, time.time() - t0), flush=True)
        del dt, P, q, Dc
        gc.collect()
    if reused:
        print('resume: %d slice files of an earlier run reused (not recomputed here)' % reused)
    print('slices s_t = %s: %d slices, %d computed here, every one a polynomial: %s; every slice computed here with '
          'positive coefficients: %s [%.0fs]' % (ka, B + 1, B + 1 - reused, good, npos_all, time.time() - t0))
    return good


def inv_vandermonde(B):
    V = flint.fmpq_mat(B + 1, B + 1, [Q(c) ** l for c in range(B + 1) for l in range(B + 1)])
    return V.inv()


def interpolate(slices, B, low, chunk=1500, log=print):
    """slices[c] = {(a, b, e): (p, q)} for c = 0..B.  Returns (N, nbad, mons): N[(l, a, b, e)] = coefficient (fmpq)
    for l <= low, nbad = number of nonzero coefficients above x1^low."""
    t0 = time.time()
    mons = set()
    for d in slices:
        mons.update(d.keys())
    mons = sorted(mons)
    Vi = inv_vandermonde(B)
    N = {}
    nbad = 0
    for c0 in range(0, len(mons), chunk):
        part = mons[c0:c0 + chunk]
        Y = flint.fmpq_mat(B + 1, len(part))
        for c, d in enumerate(slices):
            for j, mo in enumerate(part):
                v = d.get(mo)
                if v is not None:
                    Y[c, j] = Q(v[0], v[1])
        C = Vi * Y
        for j, mo in enumerate(part):
            for l in range(B + 1):
                v = C[l, j]
                if v == 0:
                    continue
                if l > low:
                    nbad += 1
                else:
                    N[(l,) + tuple(mo)] = v
        del Y, C
        gc.collect()
        log('  monomials %d..%d of %d interpolated [%.0fs]' % (c0, c0 + len(part) - 1, len(mons), time.time() - t0),
            flush=True)
    return N, nbad, mons


def mode_interp(ka):
    t0 = time.time()
    st = load_setup(ka)
    B, low = st['B'], st['L']
    slices = []
    for c in range(B + 1):
        with open(os.path.join(OUT, ka, 's_%d.pkl' % c), 'rb') as fh:
            slices.append(pickle.load(fh))
    N, nbad, mons = interpolate(slices, B, low)
    del slices
    gc.collect()
    npos = sum(1 for v in N.values() if v > 0)
    nneg = sum(1 for v in N.values() if v < 0)
    base = any(v > 0 and k[0] == 0 and k[1] == 0 and k[2] == 0 for k, v in N.items())
    base_r = any(v > 0 and k[0] == 0 and k[1] == 0 for k, v in N.items())
    degs = [max(k[i] for k in N) for i in range(4)]
    print('N~(x1, x2, r, s) for s_t = %s: %d monomials x2^a r^b s^e interpolated over x1 = 0..%d; nonzero '
          'coefficients above x1^%d: %d; %d nonzero coefficients, %d positive, %d negative; degrees (x1, x2, r, s) %s; '
          'a positive coefficient free of x1, x2, r: %s (free of x1, x2: %s) [%.0fs]'
          % (ka, len(mons), B, low, nbad, len(N), npos, nneg, degs, base, base_r, time.time() - t0))
    with open(os.path.join(OUT, 'N_%s.pkl' % ka), 'wb') as fh:
        pickle.dump({'N': {k: (int(v.p), int(v.q)) for k, v in N.items()}, 'B': B, 'L': low, 'nbad': nbad,
                     'npos': npos, 'nneg': nneg, 'base': base, 'degs': degs,
                     'ok': nbad == 0 and nneg == 0 and base}, fh, protocol=4)
    return nbad == 0 and nneg == 0 and base


def mode_final():
    import prove_gersh as PG
    t0 = time.time()
    ok = True
    for ka in ('s', '1'):
        st = load_setup(ka)
        A, OS, factors, data = setup(ka)
        same = data['facs'] == st['facs'] and data['estar'] == st['estar'] and data['bounds'] == st['bounds']
        divides = all(d <= e for d, e in zip(st['dexp'], data['estar']))
        dc_ok = all(L.one_signed_pos(f, FREE) for f, e in zip(factors, st['dexp']) if e > 0)
        B = data['bounds']['x1']
        degR = data['degDstar']['x1'] - sum(e * f.degrees()[0] for f, e in zip(factors, st['dexp']))
        low = B - degR
        consistent = (B, degR, low) == (st['B'], st['degR'], st['L'])
        print('s_t = %s: setup recomputed and equal to the saved one: %s; Den_c divides D*: %s; Den_c > 0 on the '
              'region: %s; entry denominators > 0 on the region: %s; B = %d, deg R = %d, L = %d (as saved: %s)'
              % (ka, same, divides, dc_ok, data['dens_ok'], B, degR, low, consistent), flush=True)
        ok &= same and divides and dc_ok and data['dens_ok'] and consistent
        # the first two pivots of this same matrix: equal to the conditions B_s12 / B_112 of prove_gersh.py (computed
        # there by late composition), and certified again here
        p1 = A[0][0]
        p2 = A[1][1] - A[0][1] * A[0][1] / A[0][0]
        PG.ONLY = 'B_%s12' % ka
        conds = PG.conditions('gen')
        PG.ONLY = None
        same12 = len(conds) == 2 and conds[0][1] == p1 and conds[1][1] == p2
        print('s_t = %s: pivots 1 and 2 of A equal the functions of prove_gersh.py gen only=B_%s12: %s'
              % (ka, ka, same12), flush=True)
        ok &= same12
        ok &= PG.cert('A (s_t = %s) pivot 1' % ka, p1, FREE)
        ok &= PG.cert('A (s_t = %s) pivot 2' % ka, p2, FREE)
        del A, OS, p1, p2, conds
        gc.collect()
        missing = [c for c in range(B + 1) if not os.path.exists(os.path.join(OUT, ka, 's_%d.pkl' % c))]
        with open(os.path.join(OUT, 'N_%s.pkl' % ka), 'rb') as fh:
            Nd = pickle.load(fh)
        interp_ok = Nd['ok'] and Nd['B'] == B and Nd['L'] == low and Nd['nbad'] == 0 and Nd['nneg'] == 0 and Nd['base']
        print('s_t = %s: slices missing: %s; interpolation passed (no coefficient above x1^%d, %d coefficients of N~, '
              'none negative, a positive one free of x1, x2, r): %s' % (ka, missing, low, len(Nd['N']), interp_ok),
              flush=True)
        ok &= not missing and interp_ok
    print('conclusion: det A = N~/Den_c > 0 on {x1, x2 >= 0, r > 0, s > 0} for s_t = s and s_t = 1, so the third '
          'pivot of (B_inv) is positive for all u1, u2 >= 2: %s [%.0fs]' % (ok, time.time() - t0))
    return ok


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'bounds':
        res = mode_bounds(sys.argv[2])
    elif mode == 'slices':
        res = mode_slices(sys.argv[2], resume='resume' in sys.argv[3:])
    elif mode == 'interp':
        res = mode_interp(sys.argv[2])
    elif mode == 'final':
        res = mode_final()
    else:
        raise SystemExit('unknown mode')
    print('PASS' if res else 'FAIL')
