"""The family forms of prove_gen.py against the exact numeric forms of mlib.py at members far from the sample points.
For each family given (q0,t,i with 'n' symbolic, lower bounds --imin, default 5; q0, t >= 4), the form over
Q(sizes, sigma) is computed by famlib.py as in the certificates, and N0 times its value is compared entry by entry
with mlib.reduced (exact, at rational sigma; mlib asserts its Newton window) at the members whose symbolic sizes lie
in the lists below (k up to about 80), for sigma in {3/11, 7/2, 40}.
Usage: python famcheck_mi.py FAM [FAM ...] [--imin=5]"""
import sys, time, itertools
from fractions import Fraction as Fr

opts = {'imin': '5'}
args = []
for a in sys.argv[1:]:
    if a.startswith('--'):
        k, v = a[2:].split('=')
        opts[k] = v
    else:
        args.append(a)
import famlib
famlib.SMIN['i'] = int(opts['imin'])
import mlib
from prove_fam import evalF

FAR = {'q0': (6, 17, 29), 't': (5, 19, 31), 'i': (7, 13, 19)}
T0 = time.time()
ok = True
npts = nent = 0
for a in args:
    fam = tuple(None if s == 'n' else int(s) for s in a.split(','))
    q0, t, i = fam
    C = famlib.FamCalc(q0, t, i)
    names, Q = famlib.reduced(C)
    sz = [nm for nm in C.syms if nm != 's']
    kmax = 0
    for sample in itertools.product(*[FAR[nm] for nm in sz]):
        val = dict(zip(sz, sample))
        q0v, tv, iv = val.get('q0', q0), val.get('t', t), val.get('i', i)
        if iv < famlib.SMIN['i'] and i is None:
            continue
        m = q0v + tv
        u = m + 2 * iv
        k, j = u + 1, q0v + iv
        if k > 84:
            continue
        kmax = max(kmax, k)
        N0 = famlib.norm0(C, q0v, tv, iv)
        for sg in (Fr(3, 11), Fr(7, 2), Fr(40)):
            n2, nr2, Q2 = mlib.reduced(mlib.Calc(k, j, iv, sigma=sg))
            npts += 1
            if n2 != names:
                ok = False
                print('NAMES DIFFER', a, (k, j, iv), names, n2)
                continue
            pt = list(sample) + [sg]
            for x in range(len(names)):
                for y in range(len(names)):
                    nent += 1
                    if evalF(Q[x][y], pt) * N0 != Fr(int(Q2[x][y].p), int(Q2[x][y].q)):
                        ok = False
                        print('VALUE DIFFERS', a, (k, j, iv, sg), names[x], names[y])
    print(f'{a} (i >= {famlib.SMIN["i"]} if symbolic): {len(names)} vectors, members up to k = {kmax} agree '
          f'({time.time() - T0:.0f}s)', flush=True)
print(f'{npts} (member, sigma) points, {nent} entries compared')
print('PASS' if ok else 'FAIL', f'{time.time() - T0:.0f}s')
