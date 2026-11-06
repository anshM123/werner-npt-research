"""Certificate for one block family and one two-value configuration, for every q > 0.

python cert5.py FAMILY A B [bscale] [k=5] [save=file]
  FAMILY: L3, T31, T32, T33, T42 (k5lib.py), or 'j,m'.  Configuration (A, B): the copies 0..A-1 of T0 = {0..m-1} and
  m..m+B-1 of U0 = {m..k-1} have squared weight 1, the others squared weight q.
Steps, all exact over Q(q) (ufield.RF, python-flint):
 1. the rescaled block (ylib.py, from the definitions of M, V_p, D_p, R_p and the partial transpose) over Q(q);
    every denominator of an entry is > 0 on (0, inf);
 2. (bscale = 1) the known kernel vectors (k5lib.known_kernel) are annihilated identically in q;
 3. the gauge: delete one coordinate per known kernel vector (k5lib.known_kernel);
 4. symmetric elimination of the gauge-fixed block in the natural order; every leading principal minor (the product
    of the first pivots, a reduced rational function of q) is > 0 for every q > 0 (signtest.py): numerator and
    denominator have coefficients >= 0 (not all zero); otherwise an exact evaluation at candidate points looks for a
    q where the minor is <= 0, and if there is none, an exact Sturm count decides.
Then, by Sylvester's criterion, the gauge-fixed block is positive definite for every q > 0, so the block is >= 0 with
kernel exactly the span of the known vectors (NOTES.md, Section 2).
With bscale != 1 (negative control: the B-part multiplied by bscale) step 2 is skipped and the run prints FAIL as soon
as a leading minor is not positive on (0, inf), with a rational q where it is <= 0.
Prints PASS or FAIL on the last line."""
import sys, time, pickle
from fractions import Fraction as Fr
from ufield import RF, Q, poly_positive_on_open_halfline
from signtest import positive_on_halfline, find_negative_point
import k5lib as K

t0 = time.time()


def log(*a):
    print(*a, flush=True)


def coeffs(p):
    return [str(c) for c in p.coeffs()]


def main():
    args = [a for a in sys.argv[1:] if '=' not in a]
    opts = dict(a.split('=', 1) for a in sys.argv[1:] if '=' in a)
    famname, A, B = args[0], int(args[1]), int(args[2])
    bscale = Fr(args[3]) if len(args) > 3 else Fr(1)
    k = int(opts.get('k', 5))
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
    # 1. denominators of the entries
    dens = set()
    for i in range(n):
        for l in range(i, n):
            dens.add(tuple(coeffs(M[i][l].d)))
    den_ok = True
    from ufield import P
    import flint
    for dc in dens:
        p = P([flint.fmpq(Fr(c).numerator, Fr(c).denominator) for c in dc])
        g, how = poly_positive_on_open_halfline(p)
        den_ok &= g
    log(f'  {len(dens)} distinct entry denominators, all > 0 on (0, inf): {den_ok}')
    ok &= den_ok
    # 2. kernel
    vecs, dels = K.known_kernel(k, j, T, states, s, zero, one)
    if bscale == 1:
        ann = all(all(x == 0 for x in K.matvec(M, v)) for v in vecs)
        log(f'  {len(vecs)} known kernel vectors annihilated identically in q: {ann}')
        ok &= ann
    # 3. gauge
    keep = [i for i in range(n) if i not in dels]
    G = K.principal(M, keep)
    ng = len(G)
    sym = all(G[i][l] == G[l][i] for i in range(ng) for l in range(ng))
    log(f'  gauge-fixed block: {ng} states (deleted {len(dels)}), symmetric {sym}')
    ok &= sym
    # 4. leading minors by symmetric elimination
    Aw = [list(r) for r in G]
    minors = []
    minor = one
    allpos = True
    maxn = maxd = 0
    nonneg_all = True
    failed = None
    for c in range(ng):
        p = Aw[c][c]
        if p.is_zero():
            log(f'  pivot {c + 1} vanishes identically')
            allpos = False
            failed = (c + 1, None)
            break
        minor = minor * p
        g, how, wit = positive_on_halfline(minor)
        minors.append((coeffs(minor.n), coeffs(minor.d)))
        maxn, maxd = max(maxn, minor.n.degree()), max(maxd, minor.d.degree())
        if how.count('nonnegative') < 2:
            nonneg_all = False
        log(f'     minor {c + 1}: num deg {minor.n.degree()}, den deg {minor.d.degree()}: {how}  [{time.time() - t0:.1f}s]')
        if not g:
            allpos = False
            x = wit if wit is not None else find_negative_point(minor)
            failed = (c + 1, x)
            log(f'  minor {c + 1} is not > 0 on (0, inf); nonpositive at q = {x}')
            break
        for i in range(c + 1, ng):
            if not Aw[i][c].is_zero():
                f = Aw[i][c] / p
                for l in range(c + 1, ng):
                    if not Aw[c][l].is_zero():
                        Aw[i][l] = Aw[i][l] - f * Aw[c][l]
    ok &= allpos
    log(f'  every leading principal minor > 0 on (0, inf): {allpos}; all numerators and denominators with only '
        f'nonnegative coefficients: {nonneg_all and allpos}; largest degrees: numerator {maxn}, denominator {maxd}')
    if 'save' in opts:
        with open(opts['save'], 'wb') as fh:
            pickle.dump(dict(family=famname, k=k, j=j, T=T, config=(A, B), bscale=str(bscale), states=states,
                             dels=dels, keep=keep,
                             block=[[(coeffs(x.n), coeffs(x.d)) for x in row] for row in M],
                             minors=minors), fh)
    log(f'time {time.time() - t0:.1f}s')
    log('PASS' if ok else 'FAIL')


if __name__ == '__main__':
    main()
