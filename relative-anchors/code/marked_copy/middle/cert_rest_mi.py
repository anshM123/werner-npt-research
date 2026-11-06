"""Fixed-k certificates (cert_fixed.certify: the reduced form over Q(sigma), the kernel identities, a sign-definite
gauge, every pivot of the gauge-fixed form positive on (0, inf)) for the reduced forms with KMIN <= k <= KMAX whose
whole form is not covered by a proved family (the pickles are read as in coverage_mi.py: whole-form families, or a
(u - i, i) family together with a (u - i - 1, i, 1) family).  Saves {(k, j, i): status} to OUT.
Usage: python cert_rest_mi.py KMIN KMAX OUT.pkl pickles..."""
import sys, time, pickle
from fractions import Fraction as Fr
import mlib
from cert_fixed import certify

KMIN, KMAX, OUT = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
OLDMIN = (4, 4, 2)
whole, part2, part3 = [], [], []


def norm_key(key):
    return tuple(('>=', mn) if v is None else v for v, mn in zip(key, OLDMIN))


for f in sys.argv[4:]:
    d = pickle.load(open(f, 'rb'))
    if isinstance(d, dict) and 'results' in d:
        lst = {'32': whole, '2': part2, '3': part3}[d['parts']]
        lst += [key for key, r in d['results'].items() if r.startswith('PD') or r.startswith('PART')]
    else:
        whole += [norm_key(key) for key, r in d.items() if r.startswith('PD')]


def contains(key, q0, t, i):
    for spec, v in zip(key, (q0, t, i)):
        if isinstance(spec, tuple):
            if v < spec[1]:
                return False
        elif spec != v:
            return False
    return True


T0 = time.time()
ok = True
saved = {}
for k in range(KMIN, KMAX + 1):
    t0 = time.time()
    u = k - 1
    cnt = 0
    for j in range(-1, k + 1):
        for i in range(0, u // 2 + 1):
            q0, t = j - i, u - j - i
            if any(contains(f, q0, t, i) for f in whole):
                continue
            if any(contains(f, q0, t, i) for f in part2) and any(contains(f, q0, t, i) for f in part3):
                continue
            good, info = certify(k, j, i)
            if info['status'] == 'empty':
                continue
            cnt += 1
            saved[(k, j, i)] = info['status']
            if not good:
                ok = False
                print(f'  k={k} j={j} i={i}: {info["status"]}', flush=True)
    print(f'k={k}: {cnt} reduced forms outside the proved families certified ({time.time() - t0:.1f}s)', flush=True)
    with open(OUT, 'wb') as fh:
        pickle.dump(saved, fh)
print(f'{len(saved)} forms with {KMIN} <= k <= {KMAX}; {time.time() - T0:.0f}s')
print('PASS' if ok else 'FAIL')
