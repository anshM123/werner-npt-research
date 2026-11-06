"""The symbolic forms of tlib.py (over Q(sigma)) evaluated at sigma against the exact numeric forms, and the Newton
window against all Johnson nodes.

python check_symnum.py KMIN KMAX [STRIDE]
For every k in [KMIN, KMAX] and every form of tlib.all_forms(k, 1) (every STRIDE-th form if STRIDE is given): (a) at sigma in {3/7, 5/2} every entry of the
symbolic form evaluated at sigma equals the numeric form (window 2) exactly; (b) at sigma = 3/7 the numeric form with
window 2 equals the numeric form with every genuine Johnson node (window k), exactly; (c) every denominator D of an
E^-inverse sector met in the symbolic run has nonnegative coefficients and is not zero (so it is > 0 for sigma > 0)."""
import sys, time
from fractions import Fraction as Fr
import tlib


def tofr(x):
    return Fr(int(x.p), int(x.q))


def main():
    kmin, kmax = int(sys.argv[1]), int(sys.argv[2])
    stride = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    T0 = time.time()
    ok = True
    nent = 0
    nforms = 0
    nD = 0
    for k in range(kmin, kmax + 1):
        for (j, tau, case, mode, i) in tlib.all_forms(k, 1)[::stride]:
            Cs = tlib.Calc(k, j, tau, i, case, mode, sigma=None)
            names, norms, Qs = tlib.reduced(Cs)
            if not names:
                continue
            nforms += 1
            for sk, D in Cs.D.items():
                nD += 1
                cs = D.coeffs()
                if D.is_zero() or any(c < 0 for c in cs):
                    ok = False
                    print('D NOT POSITIVE', k, j, tau, case, mode, i, sk, D, flush=True)
            for sg in (Fr(3, 7), Fr(5, 2)):
                Cn = tlib.Calc(k, j, tau, i, case, mode, sigma=sg)
                nn, _, Qn = tlib.reduced(Cn)
                assert nn == names
                for x in range(len(names)):
                    for y in range(len(names)):
                        nent += 1
                        if Qs[x][y](sg) != tofr(Qn[x][y]):
                            ok = False
                            print('SYMBOLIC != NUMERIC', k, j, tau, case, mode, i, names[x], names[y], flush=True)
            Cf = tlib.Calc(k, j, tau, i, case, mode, sigma=Fr(3, 7), window=k)
            nf, _, Qf = tlib.reduced(Cf)
            Cn = tlib.Calc(k, j, tau, i, case, mode, sigma=Fr(3, 7))
            _, _, Qn = tlib.reduced(Cn)
            for x in range(len(names)):
                for y in range(len(names)):
                    if Qf[x][y] != Qn[x][y]:
                        ok = False
                        print('WINDOW != ALL NODES', k, j, tau, case, mode, i, names[x], names[y], flush=True)
        print(f'k={k} done [{time.time() - T0:.1f}s]', flush=True)
    print(f'{nforms} forms, {nent} entries compared, {nD} sector denominators checked')
    print('PASS' if ok else 'FAIL', f'{time.time() - T0:.0f}s')


if __name__ == '__main__':
    main()
