"""Fixed-k certificates (cert_fixed.certify of marked_copy: the reduced form over Q(sigma), the kernel identities, a
sign-definite gauge, every pivot of the gauge-fixed form positive on (0, inf)) for the reduced forms with
KMIN <= k <= KMAX that no proved family contains (families read from the given result pickles, as in
coverage_mid.py).  Saves {(k, j, i): status} to OUT.
Usage: python cert_rest.py KMIN KMAX OUT.pkl results1.pkl [results2.pkl ...]"""
import sys, time, pickle, itertools
from fractions import Fraction as Fr
import mlib
from famlib import SMIN
from cert_fixed import certify

KMIN, KMAX, OUT = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
res = {}
for f in sys.argv[4:]:
    res.update(pickle.load(open(f, 'rb')))
proved = {key for key, r in res.items() if r.startswith('PD')}
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
            cands = itertools.product((q0, None) if q0 >= SMIN['q0'] else (q0,),
                                      (t, None) if t >= SMIN['t'] else (t,),
                                      (i, None) if i >= SMIN['i'] else (i,))
            if any(c in proved for c in cands):
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
