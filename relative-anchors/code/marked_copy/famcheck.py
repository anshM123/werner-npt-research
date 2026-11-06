"""The family forms of the marked-copy work (famlib.reduced, symbolic in the sizes and sigma) against the numeric forms
of mlib.py at members far from the sample points used by the certificates (which are min, min + 1, min + 3 of each
symbolic size).  Every family used by Theorems B, F and I is recomputed and evaluated at the sizes below and at two
values of sigma; N0 times the family form must equal the numeric form entry by entry.

  python famcheck.py bottom W | top W | fam1 | inv
"""
import sys, os, time, itertools
from fractions import Fraction as Fr
import flint

import famlib, mlib
from prove_fam import evalF

SIZES = {'q0': (6, 11, 23, 40), 't': (6, 11, 23, 40), 'i': (2, 6, 13)}
SIZES2 = {'q0': (5, 23), 't': (6, 31), 'i': (3, 9)}
SIGMAS = (Fr(3, 11), Fr(7, 2))


def check(fam):
    q0, t, i = fam
    if q0 is not None and t is not None and q0 + t < 0:
        return 'empty', 0
    C = famlib.FamCalc(q0, t, i)
    names, Q = famlib.reduced(C)
    if not names:
        return 'empty', 0
    sz = [nm for nm in C.syms if nm != 's']
    table = SIZES if len(sz) == 1 else SIZES2
    npts = 0
    for sample in itertools.product(*[table[nm] for nm in sz]):
        val = dict(zip(sz, sample))
        q0v, tv, iv = val.get('q0', q0), val.get('t', t), val.get('i', i)
        m = q0v + tv
        u = m + 2 * iv
        k, j = u + 1, q0v + iv
        N0 = famlib.norm0(C, q0v, tv, iv)
        for sg in SIGMAS:
            n2, nr2, Q2 = mlib.reduced(mlib.Calc(k, j, iv, sigma=sg))
            if n2 != names:
                return f'NAMES DIFFER at {(k, j, iv)}', npts
            pt = list(sample) + [sg]
            for a in range(len(names)):
                for b in range(len(names)):
                    if evalF(Q[a][b], pt) * N0 != Fr(int(Q2[a][b].p), int(Q2[a][b].q)):
                        return f'VALUES DIFFER at {(k, j, iv, sg)} {names[a]},{names[b]}', npts
            npts += 1
    return f'agrees at {npts} points', npts


def main():
    mode = sys.argv[1]
    fams = []
    if mode in ('bottom', 'top'):
        W = int(sys.argv[2])
        for b in range(0, W + 2):
            if mode == 'bottom':
                for a in range(-1, W - b + 1):
                    fams.append((a, None, b))
            else:
                for c in range(-2, W - b):
                    fams.append((None, c, b))
    elif mode == 'fam1':
        for q0 in (-1, 0, 1, 2, 3, None):
            for t in (-2, -1, 0, 1, 2, 3, None):
                for i in (0, 1, None):
                    if sum(v is None for v in (q0, t, i)) == 1:
                        fams.append((q0, t, i))
    else:
        fams.append((None, None, 0))
    T0 = time.time()
    ok = True
    tot = 0
    for fam in fams:
        t0 = time.time()
        res, npts = check(fam)
        tot += npts
        ok &= res.startswith('agrees') or res == 'empty'
        print(f'{fam}: {res} ({time.time() - t0:.1f}s)', flush=True)
    print(f'{len(fams)} families, {tot} (member, sigma) points; {time.time() - T0:.0f}s', flush=True)
    print('PASS' if ok else 'FAIL', flush=True)


if __name__ == '__main__':
    main()
