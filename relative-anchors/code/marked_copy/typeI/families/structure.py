"""The structure lemma for the family forms: dependence on (i, tau') through nu = tau' + 2i.

For each family given with i and tau' symbolic, the form over Q(q0, t, i, tau', sigma) (fixed sizes allowed for q0
and t) is rewritten with tau' = nu - 2i in Q(q0, t, i, nu, sigma): every entry must have a denominator free of i and
a numerator of degree <= DMAX in i (DMAX = 2 for the sym mode, 0 for the anti mode, where the form then depends on
(i, tau') only through nu).  Exact (python-flint fmpq_mpoly composition).
python structure.py FAM [FAM ...]"""
import sys, time
import flint
import tfamlib as T
import tprove as P


def check(case, mode, fam, C=None, names=None, Q=None):
    """(good, largest degree in i, allowed degree, number of vectors) for the family fam with i and tau' symbolic"""
    assert fam[2] is None and fam[3] is None
    if C is None:
        C = T.FamCalc(case, mode, *fam)
        names, Q = T.reduced(C)
    if True:
        old = C.syms
        new = [('nu' if nm == 'tp' else nm) for nm in old]
        ctx2 = flint.fmpq_mpoly_ctx.get(tuple(new), 'lex')
        g2 = dict(zip(new, ctx2.gens()))
        args = [g2['nu'] - 2 * g2['i'] if nm == 'tp' else g2[nm] for nm in old]
        ii = new.index('i')
        dmax = 2 if mode == 'sym' else 0
        worst = 0
        good = True
        n = len(names)
        for a in range(n):
            for b in range(n):
                x = Q[a][b]
                if x == 0:
                    continue
                num = x.n.compose(*args, ctx=ctx2)
                den = x.d.compose(*args, ctx=ctx2)
                g = num.gcd(den)
                num, den = num / g, den / g
                dn, dd = num.degrees()[ii], den.degrees()[ii]
                worst = max(worst, dn)
                if dd != 0 or dn > dmax:
                    good = False
    return good, worst, dmax, n


def dens_free_of_i(C):
    """every intermediate denominator of the calculus (C.dens), rewritten with tau' = nu - 2i, is free of i; then at a
    member (i, tau') it takes the value of the same denominator of the family with i = 2 at (2, tau' + 2i - 4)"""
    old = C.syms
    new = [('nu' if nm == 'tp' else nm) for nm in old]
    ctx2 = flint.fmpq_mpoly_ctx.get(tuple(new), 'lex')
    g2 = dict(zip(new, ctx2.gens()))
    args = [g2['nu'] - 2 * g2['i'] if nm == 'tp' else g2[nm] for nm in old]
    ii = new.index('i')
    return all(x.n.compose(*args, ctx=ctx2).degrees()[ii] == 0 for x, _ in C.dens.values()), len(C.dens)


def main():
    ok = True
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    files = [a[1:] for a in args if a.startswith('@')]
    specs = [a for a in args if not a.startswith('@')]
    for fn in files:
        specs += [s for s in open(fn).read().split() if s.split('/')[2].split(',')[2:] == ['n', 'n']]
    dens = '--dens' in sys.argv
    for s in specs:
        t0 = time.time()
        case, mode, fam = P.parse_fam(s)
        C = T.FamCalc(case, mode, *fam)
        names, Q = T.reduced(C)
        good, worst, dmax, n = check(case, mode, fam, C, names, Q)
        extra = ''
        if dens:
            gd, nd = dens_free_of_i(C)
            good &= gd
            extra = f'; its {nd} intermediate denominators free of i: {gd}'
        ok &= good
        print(f'{s}: {n} vectors; with tau\' = nu - 2i every denominator is free of i and the largest degree in i of a '
              f'numerator is {worst} (allowed {dmax}){extra}: {good} [{time.time()-t0:.1f}s]', flush=True)
    print(f'{len(specs)} forms')
    print('PASS' if ok else 'FAIL')


if __name__ == '__main__':
    main()
