"""Certificate for one block family and one two-value configuration, for every q > 0, after the isotypic reduction
of symred.py (used where the direct certificate of cert5.py is too large).

python cert_sym.py FAMILY A B [bscale] [k=6] [save=file]
Steps, exact over Q(q) (ufield.RF) and Q:
 1. the rescaled block Yr over Q(q) (ylib.py, from the definitions); every entry denominator is > 0 on (0, inf);
    den = lcm of the entry denominators, N = den Yr (polynomial entries);
 2. (bscale = 1) the known kernel vectors are annihilated identically in q;
 3. W: the G-stable complement of the known kernel K cut out by <1_y, x> = 0 for every column group y (and, for the
    letter-1 block of (R_inf)_{k-2}, <1_vac, x> = 0); a basis of W from differences and unit vectors;
 4. G = the permutations of the copies inside each of the four blocks {T on level 1, T on level q, U on level 1,
    U on level q}; the state list is G-stable and Yr is G-invariant identically in q (checked on the generating
    transpositions);
 5. for every irrep lam of G: a basis B of e_lam W (symred.py), the dimension count sum dim(lam) dim(e_lam W) = dim W,
    and F = B N B^T;
 6. every leading principal minor of every F is > 0 for every q > 0 (signtest.py: coefficient signs; otherwise a
    search for an exact witness q, then a Sturm count).
Since den > 0 on (0, inf), B Yr B^T = F/den is then positive definite for every q > 0 and every lam, so Yr|_W > 0
(symred.py) and the block is >= 0 with kernel exactly K (NOTES.md, Section 2).
With bscale != 1 (negative control) step 2 is skipped and the run prints FAIL as soon as a minor is not positive.
Prints PASS or FAIL on the last line."""
import sys, time, pickle
from fractions import Fraction as Fr
import flint
from ufield import RF, Q, P, poly_positive_on_open_halfline
from signtest import positive_on_halfline, find_negative_point
import k5lib as K
import symred as SR

t0 = time.time()


def log(*a):
    print(*a, flush=True)


def coeffs(p):
    return [str(c) for c in p.coeffs()]


def fq(x):
    x = Fr(x)
    return flint.fmpq(x.numerator, x.denominator)


def main():
    args = [a for a in sys.argv[1:] if '=' not in a]
    opts = dict(a.split('=', 1) for a in sys.argv[1:] if '=' in a)
    famname, A, B = args[0], int(args[1]), int(args[2])
    bscale = Fr(args[3]) if len(args) > 3 else Fr(1)
    k = int(opts.get('k', 6))
    j, m = K.fam(famname)
    T = tuple(range(m))
    one, zero = RF(1), RF(0)
    s = K.config_s(k, m, A, B, Q, one)
    log(f'family {famname}: (R_inf)_{j}, k = {k}, T = {T}, configuration ({A}, {B}): squared weights '
        f'{["1" if x == one else "q" for x in s]}, bscale {bscale}')
    states, M = K.build_ylib(k, j, T, s, zero, one, bscale=bscale)
    n = len(states)
    log(f'  block over Q(q): {n} states [{time.time() - t0:.1f}s]')
    ok = True
    # 1. denominators, common denominator
    den = P([1])
    dset = {}
    for i in range(n):
        for l in range(i, n):
            d = M[i][l].d
            dset[str(d)] = d
    den_ok = True
    for d in dset.values():
        g, how = poly_positive_on_open_halfline(d)
        den_ok &= g
        den = den * (d // den.gcd(d))
    den_ok &= poly_positive_on_open_halfline(den)[0]
    log(f'  {len(dset)} distinct entry denominators, all > 0 on (0, inf): {den_ok}; lcm of degree {den.degree()}')
    ok &= den_ok
    N = [[(M[i][l].n * (den // M[i][l].d)) for l in range(n)] for i in range(n)]
    # 2. kernel
    vecs, dels = K.known_kernel(k, j, T, states, s, zero, one)
    if bscale == 1:
        ann = all(all(x == 0 for x in K.matvec(M, v)) for v in vecs)
        log(f'  {len(vecs)} known kernel vectors annihilated identically in q: {ann}')
        ok &= ann
    # 3. basis of W
    groups = K.kernel_groups(k, states)
    ingroup = set(i for g in groups.values() for i in g)
    Wb = []
    for y in sorted(groups):
        g = groups[y]
        for i in g[1:]:
            Wb.append({g[0]: Fr(1), i: Fr(-1)})
    vac_k2 = []
    if m == 0 and j == k - 2:
        vac_k2 = [i for i, (a, y) in enumerate(states) if not any(a)]
        for i in vac_k2[1:]:
            Wb.append({vac_k2[0]: Fr(1), i: Fr(-1)})
    for i in range(n):
        if i not in ingroup and i not in vac_k2:
            Wb.append({i: Fr(1)})
    dimW = len(Wb)
    log(f'  W: dimension {dimW} = {n} - {len(vecs)}: {dimW == n - len(vecs)}')
    ok &= dimW == n - len(vecs)
    # 4. the group
    blocks = [list(range(0, A)), list(range(A, m)), list(range(m, m + B)), list(range(m + B, k))]
    blocks = [b_ for b_ in blocks if b_]
    mod = SR.Module(states)
    gens = []
    for b_ in blocks:
        for x in range(len(b_) - 1):
            p = list(range(k))
            p[b_[x]], p[b_[x + 1]] = b_[x + 1], b_[x]
            gens.append(tuple(p))
    inv_ok = mod.check_invariant(gens)
    for p in gens:
        im = [mod.img(p, i) for i in range(n)]
        inv_ok &= all(M[im[i]][im[l]] == M[i][l] for i in range(n) for l in range(i, n))
    log(f'  G = product of the symmetric groups of {blocks}: states G-stable and block G-invariant: {inv_ok}')
    ok &= inv_ok
    # 5. parts
    parts = []
    total = 0
    for lams in SR.irreps(blocks):
        Yg = SR.Young(k, blocks, lams)
        imgs = [mod.apply_e(Yg, v) for v in Wb]
        imgs = [v for v in imgs if v]
        Bm = SR.span_basis(imgs, n)
        d = len(Bm)
        total += Yg.dim * d
        if d == 0:
            continue
        # F = B N B^T over Q[q]
        BN = []
        for row in Bm:
            nz = [(l, fq(c)) for l, c in enumerate(row) if c != 0]
            BN.append([sum((N[i][l] * c for l, c in nz), P([0])) for i in range(n)])
        F = []
        for r1 in range(d):
            F.append([sum((BN[r1][i] * fq(c) for i, c in enumerate(Bm[r2]) if c != 0), P([0])) for r2 in range(d)])
        parts.append((lams, Yg.dim, Bm, F))
    cnt_ok = total == dimW
    log(f'  {len(parts)} parts; sum dim(lam) dim(e_lam W) = {total} = dim W: {cnt_ok} [{time.time() - t0:.1f}s]')
    ok &= cnt_ok
    # 6. minors
    saved = []
    allpos = True
    nonneg_all = True
    maxn = maxd = 0
    for lams, dl, Bm, F in parts:
        d = len(F)
        Aw = [[RF(F[r1][r2]) for r2 in range(d)] for r1 in range(d)]
        sym = all(F[r1][r2] == F[r2][r1] for r1 in range(d) for r2 in range(d))
        ok &= sym
        minor = one
        mins = []
        good = True
        for c in range(d):
            p = Aw[c][c]
            if p.is_zero():
                good = False
                log(f'     part {lams}: pivot {c + 1} vanishes identically')
                break
            minor = minor * p
            g, how, wit = positive_on_halfline(minor)
            mins.append((coeffs(minor.n), coeffs(minor.d)))
            maxn, maxd = max(maxn, minor.n.degree()), max(maxd, minor.d.degree())
            if how.count('nonnegative') < 2:
                nonneg_all = False
            if not g:
                good = False
                x = wit if wit is not None else find_negative_point(minor)
                log(f'     part {lams}: minor {c + 1} is not > 0 on (0, inf); nonpositive at q = {x}')
                break
            for i in range(c + 1, d):
                if not Aw[i][c].is_zero():
                    f = Aw[i][c] / p
                    for l in range(c + 1, d):
                        if not Aw[c][l].is_zero():
                            Aw[i][l] = Aw[i][l] - f * Aw[c][l]
        allpos &= good
        log(f'  part {lams} (dim {dl}): size {d}, symmetric {sym}, every leading minor > 0 on (0, inf): {good} '
            f'[{time.time() - t0:.1f}s]')
        saved.append(dict(lams=lams, dim=dl, B=[[str(c) for c in row] for row in Bm],
                          F=[[coeffs(x) for x in row] for row in F], minors=mins))
        if not good:
            break
    ok &= allpos
    log(f'  every leading minor of every part > 0 on (0, inf): {allpos}; all with only nonnegative coefficients: '
        f'{nonneg_all and allpos}; largest degrees: numerator {maxn}, denominator {maxd}')
    if 'save' in opts:
        with open(opts['save'], 'wb') as fh:
            pickle.dump(dict(family=famname, k=k, j=j, T=T, config=(A, B), bscale=str(bscale), states=states,
                             block=[[(coeffs(x.n), coeffs(x.d)) for x in row] for row in M],
                             den=coeffs(den), parts=saved), fh)
    log(f'time {time.time() - t0:.1f}s')
    log('PASS' if ok else 'FAIL')


if __name__ == '__main__':
    main()
