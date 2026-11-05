"""The reduced forms of l1sym.py against explicit projections of the j1lib.py block at two larger size pairs,
(u1, u2) = (5, 4) and (6, 3) (k = 9, block size 117), where every joint pattern of the pattern sums is realised with
several copies to spare, for every part (irrep on the first level; the second level by exchanging the weights).
Exact rational arithmetic, B-part unscaled.  Usage: python check_reduce_large.py"""
import time
from fractions import Fraction as Fr
from check_reduce import explicit
from l1sym import Model, reduce_part, keep_for

t0 = time.time()
ok = True
for (w0, w1, u0, u1) in [(Fr(1), Fr(2, 5), 5, 4), (Fr(3), Fr(1), 6, 3)]:
    for sw in (False, True):
        a, c = (u1, u0) if sw else (u0, u1)
        wa, wc = (w1, w0) if sw else (w0, w1)
        for part in ('inv', 'std', 'two', 'anti', 'mix'):
            if (part == 'std' and a < 2) or (part == 'two' and a < 4) or (part == 'anti' and a < 3) or \
                    (part == 'mix' and (a < 2 or c < 2)) or (part in ('inv', 'mix') and sw):
                continue
            kp, kq, g = keep_for(part, a, c)
            red = reduce_part(Model(Fr(1), [wa * wa, wc * wc], [Fr(a), Fr(c)]), part, kp, kq)
            ex = explicit(part, wa, wc, a, c, Fr(1))
            same = red['RPP'] == ex['RPP'] and red['S'] == ex['S'] and red['GP'] == ex['GP']
            if 'RPQ' in ex:
                same &= red['RPQ'] == ex['RPQ'] and red['RQQ'] == ex['RQQ']
            ok &= same
            print(f'(u1, u2) = ({u0}, {u1}), irrep on the level of size {a}: {part:4s} reduced = explicit {same} '
                  f'[{time.time() - t0:.1f}s]', flush=True)
print('PASS' if ok else 'FAIL')
