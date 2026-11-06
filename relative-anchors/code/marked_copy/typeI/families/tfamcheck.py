"""The family forms of tfamlib.py against the numeric forms of tlib.py at members far from the sanity points.

For each family given (case/mode/q0,t,i,tau' with 'n' for a symbolic size) the form is computed once over
Q(sizes, sigma); at every member whose symbolic sizes lie in the far value lists below (all combinations, k <= KMAX)
and at sigma in {3/11, 7/2, 40}, N0 times the form evaluated there must equal tlib's exact numeric form (tlib asserts
its Newton window), with the same twisted vectors.  Also checked for every family: the form is symmetric, every entry
denominator and every intermediate denominator of the calculus is one-signed after the shifts.
python tfamcheck.py KMAX FAM [FAM ...]"""
import sys, time, itertools
from fractions import Fraction as Fr
import tfamlib as T
import tprove as P

FAR = {'q0': (6, 11, 17), 't': (5, 13, 19), 'i': (3, 7, 9), 'tp': (3, 8, 12)}
SIGMAS = (Fr(3, 11), Fr(7, 2), Fr(40))


def main():
    kmax = int(sys.argv[1])
    ok = True
    npts = nent = 0
    T0 = time.time()
    for s in sys.argv[2:]:
        case, mode, fam = P.parse_fam(s)
        t0 = time.time()
        C = T.FamCalc(case, mode, *fam)
        names, Q = T.reduced(C)
        n = len(names)
        sym_ok = all(Q[a][b] == Q[b][a] for a in range(n) for b in range(n))
        den_ok = all(Q[a][b] == 0 or P.one_signed(Q[a][b].d, C) != 0 for a in range(n) for b in range(n))
        int_ok = all(P.one_signed(x.n, C) != 0 and x.d.is_constant() for x, _ in C.dens.values())
        sz = [nm for nm in C.syms if nm != 's']
        pts = []
        for sample in itertools.product(*[FAR[nm] for nm in sz]):
            vals = {nm: C.fixed[nm] for nm in T.SIZES}
            vals.update(dict(zip(sz, sample)))
            if C.member(vals)[0] <= kmax:
                pts.append((list(sample), vals))
        bad, np_ = P.sanity(C, names, Q, Fr(1), sigmas=SIGMAS, points=pts)
        good = sym_ok and den_ok and int_ok and bad is None and np_ > 0
        ok &= good
        npts += np_
        nent += np_ * n * n
        kk = [C.member(v)[0] for _, v in pts]
        print(f'{s}: {n} vectors; symmetric {sym_ok}; entry denominators one-signed {den_ok}; '
              f'{len(C.dens)} intermediate denominators one-signed {int_ok}; {np_} (member, sigma) points '
              f'with k = {min(kk) if kk else "-"}..{max(kk) if kk else "-"}: {"agree" if bad is None else bad} '
              f'[{time.time()-t0:.1f}s]', flush=True)
    print(f'{npts} (member, sigma) points, {nent} entries [{time.time()-T0:.0f}s]')
    print('PASS' if ok else 'FAIL')


if __name__ == '__main__':
    main()
