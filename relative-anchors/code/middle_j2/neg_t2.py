"""Negative controls for family 1: with the B-part scaled by 1001/1000 the gauge-fixed block (P_t deleted) is not
positive definite at explicit rational points, while the true gauge-fixed block is.

For each point: the reduced forms of the scaled gauge-fixed block (t2lev.py, invariant and standard parts) are
eliminated exactly and one of them has a nonpositive pivot; the full scaled gauge-fixed block built from the Lemma B
builder j1lib.py (bscale = 1001/1000) is not positive definite (exact LDL); the true one (bscale = 1) is.
Points: one U-level (arbitrary a, b) and two U-levels with all weights taking two values.  Usage: python neg_t2.py"""
import time
from fractions import Fraction as Fr
from j1lib import J1
from r2lib import states_T2, ldl_pd
import t2lev

t0 = time.time()
ok = True
BS = Fr(1001, 1000)


def reduced_pd(a, b, qs, us, bsc):
    alg = t2lev.point(a, b, qs, us, bsc)
    Finv, Fstd = t2lev.forms(alg, 'G')
    res = {'inv': ldl_pd(Finv)}
    for i, u in enumerate(us):
        if u >= 2:
            res[f'std{i}'] = ldl_pd(Fstd[i])
    return res


def block_pd(w, bsc):
    k = len(w)
    st = states_T2(k, 0, 1)
    M = J1(w, bscale=bsc).matrix(st)
    keep = [i for i, (a, y) in enumerate(st) if not (a is not None and a[0] == 1 and a[1] == 0)]
    return ldl_pd([[M[i][j] for j in keep] for i in keep])


# candidate points (w_t, w_tau, U-weights by level): squares of the w's are the squared weights
cands = []
for wt, wT in [(Fr(1, 30), Fr(1)), (Fr(1), Fr(1, 30)), (Fr(1, 30), Fr(1, 30)), (Fr(30), Fr(1)), (Fr(1, 10), Fr(2)),
               (Fr(1, 100), Fr(1)), (Fr(1), Fr(1))]:
    for levels in [(3, 0), (4, 0), (2, 1), (1, 2), (2, 2), (3, 1)]:
        for wq in [Fr(1, 10), Fr(10), Fr(1, 30)]:
            cands.append((wt, wT, levels, wq))
found = []
for wt, wT, (u0, u1), wq in cands:
    a, b = wt * wt, wT * wT
    qs = [Fr(1), wq * wq] if u1 else [Fr(1)]
    us = [u0, u1] if u1 else [u0]
    r_true = reduced_pd(a, b, qs, us, 1)
    r_sc = reduced_pd(a, b, qs, us, BS)
    if all(r_true.values()) and not all(r_sc.values()):
        found.append((wt, wT, (u0, u1), wq, [p for p, v in r_sc.items() if not v]))
print(f'{len(cands)} candidate points scanned with the reduced forms; {len(found)} with a failing scaled part '
      f'[{time.time() - t0:.1f}s]', flush=True)
# confirm a selection with the full block from j1lib (one level, two levels with two weight values, two levels general)
want = []
seen = set()
for f in found:
    wt, wT, (u0, u1), wq, parts = f
    key = ((u0, u1), wt == 1 and wT == 1, u1 == 0)
    if key not in seen:
        seen.add(key)
        want.append(f)
nconf = 0
for wt, wT, (u0, u1), wq, parts in want[:8]:
    w = [wt, wT] + [Fr(1)] * u0 + [wq] * u1
    pt = block_pd(w, 1)
    ps = block_pd(w, BS)
    good = pt and not ps
    nconf += good
    print(f'  w = {[str(x) for x in w]}: failing scaled parts {parts}; j1lib gauge-fixed block PD: true {pt}, '
          f'scaled {ps} -> control {"confirmed" if good else "NOT confirmed"} [{time.time() - t0:.1f}s]', flush=True)
    ok &= good
ok &= nconf >= 4
# two weight values in all (a, b in {1, q}): at least one confirmed point
tv = [f for f in found if f[0] in (Fr(1), f[3]) and f[1] in (Fr(1), f[3]) and f[2][1] > 0]
if tv:
    wt, wT, (u0, u1), wq, parts = tv[0]
    w = [wt, wT] + [Fr(1)] * u0 + [wq] * u1
    good = block_pd(w, 1) and not block_pd(w, BS)
    print(f'  two weight values: w = {[str(x) for x in w]}: failing scaled parts {parts}; confirmed with j1lib: {good}')
    ok &= good
else:
    print('  two weight values: no failing point among the candidates')
    ok = False
print('PASS' if ok else 'FAIL')
