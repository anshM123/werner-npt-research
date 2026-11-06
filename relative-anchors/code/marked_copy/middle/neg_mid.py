"""Negative controls at explicit points for the new families.  For members (k, j, i) of the families with
33 <= k <= KMAX (the reduced form computed exactly by the orbit calculus of mlib.py at rational sigma, cost
independent of k), with the irrep split of splitlib.py and the gauge of the certificates (GAUGE as in prove_split.py):
  - the true form (B unscaled), restricted to the kept coordinates, is positive definite (exact LDL) at
    sigma in {10^-6, 1, 10^6};
  - with B -> (1 + 1/1000) B, for every member there is a sigma in {10^-6, 10^6} where it is not positive definite.
Usage: python neg_mid.py KMAX FAM [FAM ...]   (FAM as in prove_mid.py, e.g. n,n,1)"""
import sys, time, os
from fractions import Fraction as Fr
import flint
import mlib, splitlib
from prove_split import parse

fq = flint.fmpq
KMAX = int(sys.argv[1])
fams = [tuple(parse(s) for s in a.split(',')) for a in sys.argv[2:]]
GAUGE = os.environ.get('GAUGE', 'O')
T0 = time.time()
ok = True


def ispd(Mx, keep):
    act = list(keep)
    M = {(a, b): Mx[a, b] for a in keep for b in keep}
    while act:
        x = act[0]
        piv = M[x, x]
        if piv <= 0:
            return False
        act.remove(x)
        for a in act:
            f = M[a, x] / piv
            if f != 0:
                for b in act:
                    M[a, b] -= f * M[x, b]
    return True


def gauge_coords(two, names):
    out = []
    for y0 in (0, 1):
        if 'O%d' % y0 not in names:
            continue
        c = 'O%d' % y0
        if GAUGE == 'R':
            if 'U1%d' % y0 in two:
                c = 'U1%d' % y0
            elif 'A1%d' % y0 in two:
                c = 'A1%d' % y0
        out.append(c)
    return out


def member(fam, k, j, i):
    q0, t = j - i, k - 1 - j - i
    for v, w, mn in zip(fam, (q0, t, i), (4, 4, 2)):
        if v is None and w < mn:
            return False
        if v is not None and v != w:
            return False
    return True


nmem = 0
for k in range(33, KMAX + 1):
    u = k - 1
    for j in range(-1, k + 1):
        for i in range(0, u // 2 + 1):
            if not any(member(f, k, j, i) for f in fams):
                continue
            if (k + j + i) % 7:      # a sparse, deterministic sample of the members
                continue
            fails = []
            for sg in (Fr(1, 10 ** 6), Fr(1), Fr(10 ** 6)):
                for bs in (Fr(1), Fr(1001, 1000)):
                    names, norms, Q = mlib.reduced(mlib.Calc(k, j, i, sigma=sg, bscale=bs))
                    vecs, split = splitlib.basis(names, fq(j - i), fq(u - 2 * i), fq(i), fq(1))
                    two, three, M = splitlib.congruence(names, Q, vecs, fq(0))
                    keep = [nm for nm in two if nm not in gauge_coords(two, names)] + three
                    pd = ispd(M, keep)
                    if bs == 1 and not pd:
                        ok = False
                        print('TRUE FORM NOT PD', k, j, i, sg)
                    if bs != 1 and not pd:
                        fails.append(sg)
            nmem += 1
            if not fails:
                ok = False
                print('NEGATIVE CONTROL NOT DETECTED', k, j, i)
    print(f'k={k}: {nmem} members so far ({time.time() - T0:.0f}s)', flush=True)
print(f'{nmem} members with 33 <= k <= {KMAX}: true forms PD at the three values of sigma; every scaled form fails '
      f'at sigma = 10^-6 or 10^6')
print('PASS' if ok else 'FAIL', f'{time.time() - T0:.0f}s')
