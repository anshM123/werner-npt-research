"""The determinant route for sym families with i and tau' symbolic (the structure lemma, as Proposition C of the README).

Family F = case/sym/[q0, t fixed or symbolic, i >= 2, tau' >= 2].  Q is the reduced form over Q(sizes, sigma) and Qg
the gauge-fixed form (one coordinate per known kernel vector deleted, as tprove.py).  Steps:
 1. preliminaries (tprove.py): Q symmetric; N0 Q equals tlib's numeric form at the sample points; every intermediate
    denominator and every entry denominator one-signed after the shifts; the kernel vectors annihilated; the gauge
    determinant a nonzero constant.  So Qg is continuous on the real region R = {sizes >= minimum, sigma > 0} and
    Q >= 0 with kernel exactly the span of the kernel vectors wherever Qg is positive definite.
 2. structure: with tau' = nu - 2i every entry of Qg is p(i, nu, ...)/d(nu, ...) with d free of i and deg_i p <= 2.
 3. D* = prod f^{e_f} over the irreducible factors f of the entry denominators (all free of i), e_f = max over
    permutations of the sum of the multiplicities (dynamic programming over subsets); T = D* det Qg is a polynomial
    with deg_i T <= B = max over permutations of sum_a deg_i p_{a, p(a)}.
 4. slices i = 2 + z, z = 0..B + 1 (the last one only for a consistency check): det Qg exactly over Q(nu, sigma)
    (q0, t fixed); D_c = prod f^{d_f} with d_f the largest multiplicity seen on the slices (d_f <= e_f is checked, so
    R = D*/D_c is a polynomial free of i); D_c det must be a polynomial P_z on every slice (exact division).
 5. interpolation in z (inverse Vandermonde over Q): N~(z, nu, sigma) of degree <= B in z with N~(z) = P_z for
    z = 0..B (with the extra node the coefficient of z^(B+1) of the interpolant must vanish).
    T - R N~ has degree <= B in i and vanishes at B + 1 values, so det Qg = N~ / D_c identically.
 6. signs: with nu = tau' + 2i, i = 2 + z, tau' = 2 + w, every coefficient of N~ is >= 0 and one free of z, w is > 0,
    and every factor of D_c is one-signed positive: det Qg > 0 on R.
 7. Qg is positive definite at one point (i, tau', sigma) = (2, 2, 1) (exact elimination).
R is connected, Qg is continuous and nonsingular on R, so Qg is positive definite on R, hence at every member and
every sigma > 0.
python tdet.py prove OUT.pkl CASE Q0,T [Q0,T ...]   (q0, t fixed; the families [q0, t, i >= 2, tau' >= 2])
python tdet.py measure CASE Q0 T               (q0, t fixed or 'n': form, structure, D*, degree bounds, one slice)"""
import sys, os, time, pickle, itertools
from fractions import Fraction as Fr
import flint
from flintfield import Field, F
import tfamlib as T
import tprove as P


def maxperm_hungarian(W):
    """the same maximum by the Hungarian method in exact integer arithmetic (O(n^3)); None if no permutation avoids
    the forbidden entries"""
    n = len(W)
    allowed = [x for row in W for x in row if x is not None]
    if not allowed:
        return None
    big = 2 * (n + 1) * (max(abs(x) for x in allowed) + 1)
    a = [[(-x if x is not None else big) for x in row] for row in W]
    INF = big * (n + 2) * 4
    u, v, p, way = [0] * (n + 1), [0] * (n + 1), [0] * (n + 1), [0] * (n + 1)
    for i in range(1, n + 1):
        p[0] = i
        j0 = 0
        minv = [INF] * (n + 1)
        used = [False] * (n + 1)
        while True:
            used[j0] = True
            i0, delta, j1 = p[j0], INF, 0
            for j in range(1, n + 1):
                if not used[j]:
                    cur = a[i0 - 1][j - 1] - u[i0] - v[j]
                    if cur < minv[j]:
                        minv[j], way[j] = cur, j0
                    if minv[j] < delta:
                        delta, j1 = minv[j], j
            for j in range(n + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minv[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while True:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
            if j0 == 0:
                break
    if any(W[p[j] - 1][j - 1] is None for j in range(1, n + 1)):
        return None
    return sum(W[p[j] - 1][j - 1] for j in range(1, n + 1))


def maxperm(W):
    """max over permutations p of sum_a W[a][p(a)] (None: forbidden), exact dynamic programming over subsets; for
    n > 18 the Hungarian method (also exact; compared with the dynamic programming in tdet.py measure)"""
    n = len(W)
    if n > 18:
        return maxperm_hungarian(W)
    best = {0: 0}
    for a in range(n):
        nb = {}
        for mask, val in best.items():
            for b in range(n):
                if mask >> b & 1 or W[a][b] is None:
                    continue
                m2 = mask | (1 << b)
                v = val + W[a][b]
                if m2 not in nb or v > nb[m2]:
                    nb[m2] = v
        best = nb
    return best.get((1 << n) - 1)


def build(case, fam, log=print):
    """the family form with the preliminaries; returns (C, names, keep, Qg) or raises AssertionError"""
    t0 = time.time()
    C = T.FamCalc(case, 'sym', *fam)
    names, Q = T.reduced(C)
    n = len(names)
    log(f'  form: {n} vectors, largest entry {max(P.tsz(Q[a][b]) for a in range(n) for b in range(n))} terms '
        f'[{time.time()-t0:.1f}s]')
    assert all(Q[a][b] == Q[b][a] for a in range(n) for b in range(n)), 'not symmetric'
    bad, npts = P.sanity(C, names, Q, Fr(1))
    assert bad is None, bad
    for key, (x, why) in C.dens.items():
        assert P.one_signed(x.n, C) != 0 and x.d.is_constant(), ('intermediate denominator', why)
    for a in range(n):
        for b in range(n):
            assert Q[a][b] == 0 or P.one_signed(Q[a][b].d, C) != 0, ('entry denominator', names[a], names[b])
    kv = T.kernel_vectors(C, names)
    for lab, v in kv:
        for a in range(n):
            s = C.zero
            for nm, cf in v.items():
                s = s + Q[a][names.index(nm)] * cf
            assert s == 0, ('kernel', lab)
    S = P.gauge_for(C, names, kv)
    assert S is not None, 'no gauge'
    keep = [c for c in range(n) if c not in S]
    Qg = [[Q[a][b] for b in keep] for a in keep]
    log(f'  preliminaries: symmetric, {npts} sanity points, {len(C.dens)} intermediate denominators and the entry '
        f'denominators one-signed, {len(kv)} kernel vectors annihilated, gauge {[names[c] for c in S]} '
        f'[{time.time()-t0:.1f}s]')
    return C, names, [names[c] for c in keep], Qg


def to_nu(C, Qg):
    """the entries with tau' = nu - 2i: ({(a, b): (num, den)}, ctx2, variable names); checks the structure"""
    old = C.syms
    new = [('nu' if nm == 'tp' else nm) for nm in old]
    ctx2 = flint.fmpq_mpoly_ctx.get(tuple(new), 'lex')
    g2 = dict(zip(new, ctx2.gens()))
    args = [g2['nu'] - 2 * g2['i'] if nm == 'tp' else g2[nm] for nm in old]
    ii = new.index('i')
    E = {}
    worst = 0
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
    assert worst <= 2, f'degree {worst} in i'
    return E, ctx2, new


def factor_data(E, n):
    facs, mult = {}, {}
    for (a, b), (num, den) in E.items():
        c, fl = den.factor()
        for f, m in fl:
            facs[str(f)] = f
            mult[(str(f), a, b)] = m
    keys = sorted(facs)
    estar = []
    for key in keys:
        W = [[(mult.get((key, a, b), 0) if (a, b) in E else None) for b in range(n)] for a in range(n)]
        estar.append(maxperm(W))
    return [facs[k] for k in keys], estar


def degree_bound(E, n, facs, estar, iv):
    W = [[((E[(a, b)][0].degrees()[iv] - E[(a, b)][1].degrees()[iv]) if (a, b) in E else None) for b in range(n)]
         for a in range(n)]
    return sum(e * f.degrees()[iv] for f, e in zip(facs, estar)) + maxperm(W)


def slice_det(E, n, ctx2, names2, fix, free):
    """det of the matrix with the variables in fix (name -> int) substituted, over Q(free)"""
    K = Field(list(free))
    gens = dict(zip(free, K.ctx.gens()))
    args = [gens[v] if v in gens else K.ctx.from_dict({(0,) * len(free): fix[v]}) for v in names2]
    A = [[(F(K, E[(a, b)][0].compose(*args, ctx=K.ctx), E[(a, b)][1].compose(*args, ctx=K.ctx))
           if (a, b) in E else K(0)) for b in range(n)] for a in range(n)]
    det = K(1)
    for c in range(n):
        if A[c][c] == 0:
            r = next((r for r in range(c + 1, n) if A[r][r] != 0), None)
            if r is None:
                r = next(r for r in range(c + 1, n) if A[r][c] != 0)
                A[c], A[r] = A[r], A[c]
                det = -det
            else:
                A[c], A[r] = A[r], A[c]
                for row in A:
                    row[c], row[r] = row[r], row[c]
        p = A[c][c]
        det = det * p
        for r in range(c + 1, n):
            if A[r][c] != 0:
                f = A[r][c] / p
                for cc in range(c + 1, n):
                    if A[c][cc] != 0:
                        A[r][cc] = A[r][cc] - f * A[c][cc]
    return K, det


def multiplicity(f, Pl):
    e = 0
    while True:
        q, r = divmod(Pl, f)
        if r != 0:
            return e
        Pl = q
        e += 1


def prove(case, q0, t, out=None, log=print):
    T0 = time.time()
    fam = (q0, t, None, None)
    C, names, keep, Qg = build(case, fam, log)
    n = len(keep)
    E, ctx2, names2 = to_nu(C, Qg)
    log(f'  structure: every entry of the gauge-fixed form ({n} x {n}) is p/d with d free of i and deg_i p <= 2 '
        f'[{time.time()-T0:.1f}s]')
    facs, estar = factor_data(E, n)
    ii = names2.index('i')
    W = [[(E[(a, b)][0].degrees()[ii] if (a, b) in E else None) for b in range(n)] for a in range(n)]
    B = maxperm(W)
    log(f'  {len(facs)} irreducible factors of the entry denominators, exponents of D* {estar}; '
        f'deg_i T <= B = {B} [{time.time()-T0:.1f}s]')
    free = [v for v in names2 if v not in ('i',)]
    assert free == ['nu', 's'], free
    slices = []
    obs = [0] * len(facs)
    for z in range(B + 2):
        t1 = time.time()
        K, det = slice_det(E, n, ctx2, names2, {'i': 2 + z}, free)
        assert det != 0, 'zero determinant on a slice'
        for k_, f in enumerate(facs):
            fs = f.compose(*[K.ctx.from_dict({(0, 0): 2 + z}) if v == 'i' else K.ctx.gens()[free.index(v)]
                             for v in names2], ctx=K.ctx)
            m_ = multiplicity(fs, det.d)
            obs[k_] = max(obs[k_], m_)
        slices.append(det)
        if z % 5 == 0:
            log(f'    slice i = {2 + z}: det {len(det.n)}/{len(det.d)} terms [{time.time()-t1:.1f}s]')
    assert all(o <= e for o, e in zip(obs, estar)), ('D_c does not divide D*', obs, estar)
    # D_c and the slice polynomials
    K = Field(free)
    Pz = []
    for z, det in enumerate(slices):
        Dc = K.ctx.from_dict({(0, 0): 1})
        for f, o in zip(facs, obs):
            if o:
                fs = f.compose(*[K.ctx.from_dict({(0, 0): 2 + z}) if v == 'i' else K.ctx.gens()[free.index(v)]
                                 for v in names2], ctx=K.ctx)
                Dc = Dc * fs ** o
        qq, r = divmod(Dc, det.d)
        assert r == 0, f'D_c det not a polynomial on slice {z}'
        Pz.append(det.n * qq)
    log(f'  D_c exponents {obs}; every slice gives a polynomial D_c det [{time.time()-T0:.1f}s]')
    # interpolation in z
    mons = sorted(set().union(*[set(p.to_dict().keys()) for p in Pz]))
    # B + 2 nodes: the interpolant has degree <= B + 1 and its coefficient of z^(B+1) must vanish (a consistency check
    # of the bound B; the proof needs only deg_i T <= B and B + 1 of the nodes)
    V = flint.fmpq_mat(B + 2, B + 2, [flint.fmpq(z) ** l for z in range(B + 2) for l in range(B + 2)]).inv()
    Y = flint.fmpq_mat(B + 2, len(mons))
    dicts = [p.to_dict() for p in Pz]
    for z, d in enumerate(dicts):
        for jx, mo in enumerate(mons):
            c = d.get(mo)
            if c is not None:
                Y[z, jx] = c
    Cf = V * Y
    ctx3 = flint.fmpq_mpoly_ctx.get(('z', 'nu', 's'), 'lex')
    terms = {}
    top = sum(1 for jx in range(len(mons)) if Cf[B + 1, jx] != 0)
    assert top == 0, f'{top} nonzero coefficients of z^(B+1): the degree bound is violated'
    for jx, mo in enumerate(mons):
        for l in range(B + 1):
            c = Cf[l, jx]
            if c != 0:
                terms[(l,) + tuple(mo)] = c
    Nt = ctx3.from_dict(terms)
    log(f'  interpolation over z = 0..{B + 1}: the coefficient of z^{B + 1} vanishes for all {len(mons)} monomials')
    log(f'  N~(z, nu, s): {len(Nt)} terms, degrees {Nt.degrees()} [{time.time()-T0:.1f}s]')
    # signs: nu = tau' + 2i = (2 + w) + 2(2 + z)
    ctx4 = flint.fmpq_mpoly_ctx.get(('z', 'w', 's'), 'lex')
    zz, ww, ss = ctx4.gens()
    Ns = Nt.compose(zz, 6 + 2 * zz + ww, ss, ctx=ctx4)
    co = Ns.to_dict()
    nneg = sum(1 for c in co.values() if c < 0)
    base = any(c > 0 and e[0] == 0 and e[1] == 0 for e, c in co.items())
    log(f'  N~ with nu = tau\' + 2i, i = 2 + z, tau\' = 2 + w: {len(co)} terms, {nneg} negative coefficients, a '
        f'positive one free of z, w: {base} [{time.time()-T0:.1f}s]')
    fpos = []
    for f, o in zip(facs, obs):
        if o:
            g = f.compose(*[2 + zz if v == 'i' else (6 + 2 * zz + ww if v == 'nu' else ss) for v in names2], ctx=ctx4)
            gd = g.to_dict()
            pos = all(c > 0 for c in gd.values()) and any(e[0] == 0 and e[1] == 0 for e in gd)
            fpos.append(pos)
    log(f'  every factor of D_c positive on the region: {all(fpos)}')
    # positive definite at one point
    pt = {'i': 2, 'tp': 2, 's': 1}
    vals = [flint.fmpq(pt[v]) for v in C.syms]
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
    log(f'  gauge-fixed form positive definite at (i, tau\', sigma) = (2, 2, 1): {pd} [{time.time()-T0:.1f}s]')
    ok = nneg == 0 and base and all(fpos) and pd
    info = {'names': names, 'keep': keep, 'B': B, 'estar': estar, 'dexp': obs, 'Nterms': len(co),
            'Ndeg': Nt.degrees(), 'time': time.time() - T0}
    if out:
        with open(out, 'wb') as fh:
            pickle.dump({'ok': ok, 'info': info, 'N': {e: (int(c.p), int(c.q)) for e, c in co.items()}}, fh)
    return ok, info


def measure(case, fam, log=print):
    T0 = time.time()
    C, names, keep, Qg = build(case, fam, log)
    n = len(keep)
    E, ctx2, names2 = to_nu(C, Qg)
    big = max(len(nm) for nm, dn in E.values())
    bigd = max(len(dn) for nm, dn in E.values())
    log(f'  structure holds; entries in {names2}: numerators up to {big} terms, denominators up to {bigd} terms '
        f'[{time.time()-T0:.1f}s]')
    facs, estar = factor_data(E, n)
    log(f'  {len(facs)} irreducible factors of the entry denominators; exponents of D* {estar} [{time.time()-T0:.1f}s]')
    bounds = {}
    for iv, v in enumerate(names2):
        if v == 'i':
            W = [[(E[(a, b)][0].degrees()[iv] if (a, b) in E else None) for b in range(n)] for a in range(n)]
            bounds[v] = maxperm(W)
            if n <= 18:
                assert maxperm_hungarian(W) == bounds[v], "Hungarian method and dynamic programming disagree"
        else:
            bounds[v] = degree_bound(E, n, facs, estar, iv)
            if n <= 18:
                Wv = [[((E[(a, b)][0].degrees()[iv] - E[(a, b)][1].degrees()[iv]) if (a, b) in E else None)
                       for b in range(n)] for a in range(n)]
                assert maxperm_hungarian(Wv) == maxperm(Wv), "Hungarian method and dynamic programming disagree"
    log(f'  degree bounds of T = D* det Qg: {bounds} [{time.time()-T0:.1f}s]')
    fix = {v: (C.smin['i'] if v == 'i' else C.smin[v]) for v in names2 if v in ('q0', 't', 'i')}
    free = ['nu', 's']
    t1 = time.time()
    K, det = slice_det(E, n, ctx2, names2, fix, free)
    log(f'  one slice ({fix} fixed, det over Q(nu, s)): {len(det.n)}/{len(det.d)} terms, {time.time()-t1:.1f}s')
    return bounds, time.time() - t1


def measure_plain(case, fam, log=print):
    """the same measurements without the structure lemma (i or tau' fixed): degree bounds of T = D* det Qg in each
    size variable and the time of one slice with every symbolic size but the first fixed at its minimum"""
    T0 = time.time()
    C, names, keep, Qg = build(case, fam, log)
    n = len(keep)
    E = {(a, b): (x.n, x.d) for a, row in enumerate(Qg) for b, x in enumerate(row) if x != 0}
    names2 = C.syms
    facs, estar = factor_data(E, n)
    log(f'  {len(facs)} irreducible factors of the entry denominators; exponents of D* {estar} [{time.time()-T0:.1f}s]')
    bounds = {v: degree_bound(E, n, facs, estar, iv) for iv, v in enumerate(names2)}
    log(f'  degree bounds of T = D* det Qg: {bounds} [{time.time()-T0:.1f}s]')
    szs = [v for v in names2 if v != 's']
    free = [szs[0], 's']
    fix = {v: C.smin[v] for v in szs[1:]}
    t1 = time.time()
    K, det = slice_det(E, n, None, names2, fix, free)
    log(f'  one slice ({fix} fixed, det over Q({szs[0]}, s)): {len(det.n)}/{len(det.d)} terms, {time.time()-t1:.1f}s')
    return bounds, time.time() - t1


if __name__ == '__main__':
    mode = sys.argv[1]
    logf = lambda m: print(m, flush=True)
    if mode == 'measure_plain':
        case, mode2, fam = P.parse_fam(sys.argv[2])
        print(P.fam_label(case, mode2, fam), flush=True)
        measure_plain(case, fam, logf)
        print('PASS')
        sys.exit(0)
    if mode == 'prove':
        out, case = sys.argv[2], sys.argv[3]
        T0 = time.time()
        results, infos = {}, {}
        allok = True
        for spec in sys.argv[4:]:
            q0, t = (int(x) for x in spec.split(','))
            lab = P.fam_label(case, 'sym', (q0, t, None, None))
            print(lab, flush=True)
            key = (case, 'sym', q0, t, ('>=', T.SMIN['i']), ('>=', T.SMIN['tp']))
            try:
                ok, info = prove(case, q0, t, None, logf)
                res = 'PD on complement (determinant route)' if ok else 'FAILED (a sign or the base point)'
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
    else:
        case = sys.argv[2]
        q0 = None if sys.argv[3] == 'n' else int(sys.argv[3])
        t = None if sys.argv[4] == 'n' else int(sys.argv[4])
        print(P.fam_label(case, 'sym', (q0, t, None, None)), flush=True)
        measure(case, (q0, t, None, None), logf)
        print('PASS')
