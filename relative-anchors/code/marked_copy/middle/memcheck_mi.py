"""Exact checks at explicit members (k, j, i) with KMIN <= k <= KMAX of the families given (a sparse deterministic
sample: k + j + i divisible by STEP), from the numeric orbit calculus of mlib.py (exact rationals, any k):
  - the known kernel vectors (mlib.kernel_vectors) are annihilated by Q_i;
  - with the split of splitlib.py the off-diagonal block vanishes, and the gauge-fixed form (the coordinates
    U_{1,y0} / A_{1,y0} deleted, as in the certificates) is positive definite by exact LDL, for sigma in
    {10^-6, 1/7, 1, 3, 10^6}: so Q_i >= 0 with kernel exactly the span of the known kernel vectors at these points;
    the three-row part (the S vectors) alone is positive definite there as well;
  - negative control, B -> (1 + 1/1000) B: the gauge-fixed form is not positive definite at sigma = 10^-6 or 10^6;
    this is required at every sampled member of the generic family [q0 >= 4, t >= 4, i >= IMIN] and counted per
    family for the others (where the gap can exceed the change of B); the members where the three-row part alone
    fails are counted as well.
A family is written q0,t,i with 'n' symbolic; the lower bounds are q0 >= 4, t >= 4 and i >= IMIN (--imin, default 5).
Usage: python memcheck_mi.py KMIN KMAX STEP FAM [FAM ...] [--imin=5]"""
import sys, time
from fractions import Fraction as Fr
import flint
import mlib, splitlib

fq = flint.fmpq
opts = {'imin': '5'}
args = []
for a in sys.argv[1:]:
    if a.startswith('--'):
        k, v = a[2:].split('=')
        opts[k] = v
    else:
        args.append(a)
KMIN, KMAX, STEP = int(args[0]), int(args[1]), int(args[2])
IMIN = int(opts['imin'])
fams = [tuple(None if s == 'n' else int(s) for s in a.split(',')) for a in args[3:]]
T0 = time.time()
ok = True


def ispd(M, keep):
    act = list(keep)
    M = {(a, b): M[a, b] for a in keep for b in keep}
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


def member(fam, k, j, i):
    q0, t = j - i, k - 1 - j - i
    for v, w, mn in zip(fam, (q0, t, i), (4, 4, IMIN)):
        if v is None and w < mn:
            return False
        if v is not None and v != w:
            return False
    return True


def gauge_coords(two, names):
    out = []
    for y0 in (0, 1):
        if 'O%d' % y0 not in names and 'A0%d' % y0 not in names and 'A1%d' % y0 not in names:
            continue
        for c in ('U1%d' % y0, 'A1%d' % y0, 'O%d' % y0):
            if c in two:
                out.append(c)
                break
    return out


nmem = nneg = n3fail = 0
byfam = {}
for k in range(KMIN, KMAX + 1):
    u = k - 1
    for j in range(-1, k + 1):
        for i in range(0, u // 2 + 1):
            if not any(member(f, k, j, i) for f in fams) or (k + j + i) % STEP:
                continue
            fails, fails3 = [], []
            for sg in (Fr(1, 10 ** 6), Fr(1, 7), Fr(1), Fr(3), Fr(10 ** 6)):
                for bs in (Fr(1), Fr(1001, 1000)):
                    if bs != 1 and sg not in (Fr(1, 10 ** 6), Fr(10 ** 6)):
                        continue
                    names, norms, Q = mlib.reduced(mlib.Calc(k, j, i, sigma=sg, bscale=bs))
                    vecs, split = splitlib.basis(names, fq(j - i), fq(u - 2 * i), fq(i), fq(1))
                    two, three, M = splitlib.congruence(names, Q, vecs, fq(0))
                    if any(M[a, b] != 0 for a in two for b in three):
                        ok = False
                        print('SPLIT NOT BLOCK DIAGONAL', k, j, i, sg, bs)
                    if bs == 1:
                        sgq = fq(sg.numerator, sg.denominator)
                        kv = mlib.kernel_vectors(names, k, j, i, sgq)
                        for lab, v in kv:
                            for a in range(len(names)):
                                s = fq(0)
                                for nm, cf in v.items():
                                    s += Q[a][names.index(nm)] * cf
                                if s != 0:
                                    ok = False
                                    print('KERNEL VECTOR NOT ANNIHILATED', k, j, i, sg, lab)
                        g = gauge_coords(two, names)
                        if len(g) != len(kv):
                            ok = False
                            print('GAUGE DOES NOT MATCH THE KERNEL', k, j, i, g, [lab for lab, _ in kv])
                    keep = [nm for nm in two if nm not in gauge_coords(two, names)] + three
                    pd = ispd(M, keep)
                    pd3 = ispd(M, three)
                    if bs == 1 and not (pd and pd3):
                        ok = False
                        print('TRUE FORM NOT POSITIVE DEFINITE ON THE COMPLEMENT', k, j, i, sg)
                    if bs != 1:
                        if not pd:
                            fails.append(sg)
                        if not pd3:
                            fails3.append(sg)
            nmem += 1
            fam = [f for f in fams if member(f, k, j, i)][0]
            byfam[fam] = byfam.get(fam, [0, 0, 0])
            byfam[fam][0] += 1
            byfam[fam][1] += bool(fails)
            byfam[fam][2] += bool(fails3)
            if fails:
                nneg += 1
            elif fam == (None, None, None):
                ok = False
                print('NEGATIVE CONTROL NOT DETECTED (generic family)', k, j, i)
            else:
                print('negative control not detected at an edge member', k, j, i, (j - i, k - 1 - j - i, i))
            if fails3:
                n3fail += 1
    print(f'k={k}: {nmem} members so far ({time.time() - T0:.0f}s)', flush=True)
lab = lambda f: '[' + ', '.join(f'{nm}>={mn}' if v is None else f'{nm}={v}'
                                for v, nm, mn in zip(f, ('q0', 't', 'i'), (4, 4, IMIN))) + ']'
for f, (a, b, c) in byfam.items():
    print(f'  {lab(f)}: {a} members; with B scaled the gauge-fixed form fails at one of sigma = 10^-6, 10^6 for {b}, '
          f'the three-row part alone for {c}')
print(f'{nmem} members with {KMIN} <= k <= {KMAX}: Q_i >= 0 with kernel exactly the known one at five values of '
      f'sigma (exact); with B scaled the gauge-fixed form fails at sigma = 10^-6 or 10^6 for {nneg} of them (required '
      f'for every member of the generic family), its three-row part alone for {n3fail}')
print('PASS' if ok else 'FAIL', f'{time.time() - T0:.0f}s')
