"""Exploration: locate rational witnesses for the negative controls (no claim rests on it)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
import prove_em as P
Q = flint.fmpq
def val(rf, pt):
    pt = [Q(v) for v in pt]
    return rf.n(*pt) / rf.d(*pt)
for H in (Q(1, 4), Q(4, 5)):
    _, _, conds = P.run(H)
    for name, pts in (('A_inv', [(1, 1, 0, Q(1, 10**4), 0), (1, 1, 0, Q(1, 100), 0), (1, 1, 0, Q(1, 10**4), 10)]),
                      ('B_inv_det', [(Q(1, 1000), 1, 0, Q(1, 10**4), 0), (Q(1, 100), 0, 1, Q(1, 100), 0), (Q(1, 1000), 0, 1, Q(1, 10**4), 0)]),
                      ('B_std', [(Q(1, 1000), 0, 1, Q(1, 10**4), 0), (Q(1, 1000), 0, 1, 100, 0), (Q(1, 100), 0, 1, 1, 0)]),
                      ('B_inv_00', [(Q(1, 1000), 0, 1, Q(1, 10**4), 0), (Q(1, 1000), 0, 1, 100, 0)])):
        for pt in pts:
            print(H, name, pt, float(val(conds[name], pt)))
