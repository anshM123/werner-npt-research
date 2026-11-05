"""Negative controls at explicit points, exact.  With the B-part scaled, B -> (1 + 1/1000) B (bscale of j1lib.py), the
gauge-fixed block (the states (min X, X) deleted) built by j1lib.py is not positive definite, while the true block is
(exact LDL), at an explicit two-level rational point for every two-level size pair with k = 5, 6, 7 (u0 >= u1 >= 1).
At the same point the reduced forms (fixedred.py with bscale 1001/1000) name the parts with a non-positive pivot.
Usage: python neg_points.py"""
import time
from fractions import Fraction as Fr
from j1lib import J1, mu
from r2lib import states_L1, ldl_pd
from fixedred import build_block, reduce_all, leading_pivots

t0 = time.time()
ok = True
CAND = [Fr(1, 30), Fr(10), Fr(30), Fr(1, 10)]


def gauge_pd(u0, u1, w1, b):
    k = u0 + u1
    w = [Fr(1)] * u0 + [w1] * u1
    st = states_L1(k)
    M = J1(w, bscale=b).matrix(st)
    mus = [mu(x, w) for x in st]
    keep = [i for i, (a, y) in enumerate(st) if not (a is not None and a[1] == min(c for c in range(k) if y[c] == 1))]
    return ldl_pd([[M[i][j] * mus[i] * mus[j] for j in keep] for i in keep])


for k in (5, 6, 7):
    for u1 in range(1, k // 2 + 1):
        u0 = k - u1
        found = None
        for w1 in CAND:
            true_pd = gauge_pd(u0, u1, w1, Fr(1))
            scaled_pd = gauge_pd(u0, u1, w1, Fr(1001, 1000))
            if true_pd and not scaled_pd:
                found = w1
                break
        if found is None:
            ok = False
            print(f'k={k} (u0, u1) = ({u0}, {u1}): no failing point among {[str(x) for x in CAND]}', flush=True)
            continue
        q0 = found * found
        parts, _ = reduce_all(build_block(u0, u1, q0, Fr(1), Fr(1001, 1000)), u0, u1, Fr(1001, 1000))
        bad = []
        for r in parts:
            piv = leading_pivots(r['S'], Fr(1)) if r['S'] else []
            if any(p is None or p <= 0 for p in piv):
                bad.append(r['lam'])
        tparts, _ = reduce_all(build_block(u0, u1, q0, Fr(1)), u0, u1)
        tgood = all(all(p is not None and p > 0 for p in leading_pivots(r['S'], Fr(1))) for r in tparts if r['S'])
        ok &= bool(bad) and tgood
        print(f'k={k} (u0, u1) = ({u0}, {u1}), w = {found} (q = {q0}): gauge-fixed block PD for the true B, not PD for '
              f'(1 + 1/1000)B; reduced forms: true all PD {tgood}, scaled parts with a non-positive pivot: {bad} '
              f'[{time.time() - t0:.1f}s]', flush=True)
print('PASS' if ok else 'FAIL')
