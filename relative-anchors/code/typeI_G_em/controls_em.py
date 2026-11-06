"""Negative controls for Theorem EM (exact).
 1. The certificate of prove_em.py fails (negative coefficients) with the constant moved to H = 1/4 and to H = 4/5.
 2. Explicit rational witnesses, computed from the exact G of the whole block (Lemma B builder, no closed form), that
    the inequalities themselves are false there: (A_inv) with H = 1/4 at w_T = (1, 1), om = 1/10, u = 3; (B_std) with
    H = 4/5 at w_T = (1/10, 1), om = 1, u = 3; (B_inv_det) with H = 4/5 at w_T = (1/10, 1), om = 1/10, u = 3.
 3. B -> (1 + 1/1000) B: at explicit rational configurations with equal U-weights the gauge-deleted block (P_{r0,h}
    removed) has a negative LDL^T pivot, while the true block (B unscaled) has all pivots positive there; and the
    Lemma B builder with the scaled B differs from the closed form (so the comparison of check_e2e.py is not vacuous).
python controls_em.py"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from emlib import direct_G, direct_conditions, Q
from cf import TypeI
import prove_em as P

ok = True
t0 = time.time()
for H in (Q(1, 4), Q(4, 5)):
    allok, res, _ = P.run(H)
    failed = [k for k, v in res.items() if not v]
    print('certificate with H = %s: certified %s; failing conditions %s' % (H, allok, failed), flush=True)
    ok &= not allok
wit = [(Q(1, 4), 'A_inv', (Fr(1), Fr(1)), Fr(1, 10), 3),
       (Q(4, 5), 'B_std', (Fr(1, 10), Fr(1)), Fr(1), 3),
       (Q(4, 5), 'B_inv_det', (Fr(1, 10), Fr(1)), Fr(1, 10), 3)]
for H, name, wT, om, u in wit:
    dc = direct_conditions(wT, om, u, H)
    vals = dc[name]
    neg = min(vals)
    print('witness H = %s, %s at w_T = %s, om = %s, u = %d: values %s -> %s' % (
        H, name, [str(x) for x in wT], om, u, [float(v) for v in vals], 'FALSE there (as it must be)' if neg < 0 else
        'not a witness'), flush=True)
    ok &= neg < 0
    dc2 = direct_conditions(wT, om, u, Q(1, 2))
    ok &= all(v > 0 for vals2 in dc2.values() for v in vals2)
found = 0
for wT, om, u in (((Fr(1, 100), Fr(1, 100)), Fr(1), 3), ((Fr(1), Fr(1, 100)), Fr(1, 100), 3),
                  ((Fr(1, 30), Fr(1, 20), Fr(1, 10)), Fr(1), 4)):
    _, piv1, M1, lab = direct_G(wT, om, u)
    _, pivs, Ms, _ = direct_G(wT, om, u, bscale=Fr(1001, 1000))
    m = len(wT)
    w = list(wT) + [om] * u
    Mc = TypeI(w, list(range(m))).build()
    n = len(lab)
    differs = any(Ms[i][j] != Mc[i][j] for i in range(n) for j in range(n))
    negs = sum(1 for p in pivs if p < 0)
    print('B -> (1001/1000) B at w_T = %s, om = %s, u = %d: negative pivots %d (true block: %d negative, min pivot %.3e); '
          'scaled builder differs from closed form: %s' % ([str(x) for x in wT], om, u, negs,
                                                            sum(1 for p in piv1 if p <= 0), float(min(piv1)), differs),
          flush=True)
    ok &= all(p > 0 for p in piv1) and differs
    found += negs > 0
ok &= found > 0
print('time %.1f s' % (time.time() - t0))
print('PASS' if ok else 'FAIL')
