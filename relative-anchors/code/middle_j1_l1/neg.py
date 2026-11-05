"""Negative controls: with the B-part scaled, B -> (1 + 1/1000) B (j1lib bscale; in l1sym.Model the parameter b),
the gauge-fixed letter-1 block (diagonal states P_{c,c} deleted, the gauge of the certificates) is not positive
definite at explicit two-level rational weight vectors, while the true one is.  At each point the reduced forms of
l1sym.py give the same verdicts (for b = 1 and for b = 1001/1000), and the parts that fail are listed.  Among the
points are three of the family u1 = 1 (one copy with its own weight).  The symbolic certificates with the scaled
B-part are run (and must fail) from run_all.sh.
Usage: python neg.py"""
import time
from fractions import Fraction as Fr
from check_block import block_verdict, parts_at, ldl_pd
from l1sym import Model, reduce_part, keep_for

t0 = time.time()
b1 = Fr(1001, 1000)


def failing_parts(w0, w1, u0, u1, b):
    out = []
    for part, sw in parts_at(u0, u1):
        a, c = (u1, u0) if sw else (u0, u1)
        qa, qc = (w1 * w1, w0 * w0) if sw else (w0 * w0, w1 * w1)
        kp, kq, gauge = keep_for(part, a, c)
        red = reduce_part(Model(Fr(1), [qa, qc], [Fr(a), Fr(c)], b), part, kp, kq)
        S = red['S']
        keep = [i for i in range(len(S)) if i not in gauge]
        good = ldl_pd([[S[i][j] for j in keep] for i in keep]) and (ldl_pd(red['RQQ']) if 'RQQ' in red else True)
        if not good:
            out.append(part + ('(level 1)' if sw else ''))
    return out


ok = True
pts = [(2, 2, Fr(1, 10)), (3, 1, Fr(1, 10)), (1, 3, Fr(1, 10)), (1, 4, Fr(10)), (4, 1, Fr(1, 10)), (3, 2, Fr(10)),
       (2, 3, Fr(1, 30))]
for (u0, u1, w1) in pts:
    w = [Fr(1)] * u0 + [w1] * u1
    k0, pd_true = block_verdict(w)
    k1, pd_scaled = block_verdict(w, b1)
    f_true = failing_parts(Fr(1), w1, u0, u1, Fr(1))
    f_scaled = failing_parts(Fr(1), w1, u0, u1, b1)
    agree = (pd_true == (not f_true)) and (pd_scaled == (not f_scaled))
    good = pd_true and not pd_scaled and agree
    ok &= good
    print(f'(u0, u1) = ({u0}, {u1}), w = 1 on level 0, {w1} on level 1: true gauge-fixed block PD {pd_true}; '
          f'scaled gauge-fixed block PD {pd_scaled}; reduced forms agree {agree}; failing parts (scaled): {f_scaled} '
          f'[{time.time() - t0:.1f}s]', flush=True)
# the failure depends on the gauge: in the gauge of middle_j1/check_l1.py (the states with row r = 0 deleted) the
# scaled block is still positive definite at w = (1, 1, 1/10, 1/10) and fails at (1, 1, 1/30, 1/30), the point
# reported there; in the diagonal gauge used here it fails at both
from fractions import Fraction as Fr
from j1lib import J1, block_states, mu


def gauge_row0_pd(w, b):
    k = len(w)
    J = J1(w, bscale=b)
    st = block_states(k, ())
    M = J.matrix(st)
    m = [mu(z, J.w) for z in st]
    Mt = [[m[i] * M[i][j] * m[j] for j in range(len(st))] for i in range(len(st))]
    keep = [i for i, (a, y) in enumerate(st) if not (a is not None and a[1] == 0)]
    return ldl_pd([[Mt[i][j] for j in keep] for i in keep])


for d in (10, 30):
    w = [Fr(1), Fr(1), Fr(1, d), Fr(1, d)]
    r0_true, r0_scaled = gauge_row0_pd(w, Fr(1)), gauge_row0_pd(w, b1)
    dg_scaled = block_verdict(w, b1)[1]
    print(f'w = (1, 1, 1/{d}, 1/{d}): row-0 gauge: true PD {r0_true}, scaled PD {r0_scaled}; diagonal gauge: scaled PD '
          f'{dg_scaled}', flush=True)
    ok &= r0_true and (r0_scaled == (d == 10)) and not dg_scaled
print('PASS' if ok else 'FAIL')
