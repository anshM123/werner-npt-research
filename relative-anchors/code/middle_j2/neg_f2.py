"""Negative controls for family 2: with the B-part scaled by 1001/1000 the gauge-fixed block (the diagonal states
P_{c,c} deleted) is not positive definite at explicit rational points, while the true gauge-fixed block is.

For each point: the reduced forms of the scaled block (f2sym.py, every part, gauge coordinates deleted) are eliminated
exactly and one of them has a nonpositive pivot; the full scaled gauge-fixed block built from the Lemma B builder
j1lib.py (bscale = 1001/1000) is not positive definite (exact LDL), and the true one (bscale = 1) is.
Usage: python neg_f2.py"""
import time
from fractions import Fraction as Fr
from j1lib import J1
from r2lib import states_T1, ldl_pd
from f2sym import Model2, reduce_part2, assemble

t0 = time.time()
ok = True
BS = Fr(1001, 1000)


def parts_at(u0, u1):
    out = [('inv', False)]
    for sw in (False, True):
        a, b = (u1, u0) if sw else (u0, u1)
        if a >= 2:
            out.append(('std', sw))
        if a >= 4:
            out.append(('two', sw))
        if a >= 3:
            out.append(('anti', sw))
    if u0 >= 2 and u1 >= 2:
        out.append(('mix', False))
    return out


def reduced_pd(a, q0, q1, u0, u1, bsc):
    res = {}
    for part, sw in parts_at(u0, u1):
        Q = [q1, q0] if sw else [q0, q1]
        U = [u1, u0] if sw else [u0, u1]
        M = Model2(Fr(1), a, Q, U, bsc)
        red = reduce_part2(M, part, min(U[0], 3), min(U[1], 3))
        F, iR, iP, iQ = assemble(red)
        keep = [i for i in range(len(F)) if i not in [iP[g] for g in red['gauge']]]
        res[part + ('[1]' if sw else '')] = ldl_pd([[F[i][j] for j in keep] for i in keep])
    return res


def block_pd(w, bsc):
    k = len(w)
    st = states_T1(k, 0)
    M = J1(w, bscale=bsc).matrix(st)
    keep = [i for i, (a, y) in enumerate(st) if not (a is not None and a[0] == 1 and y[a[1]] == 1)]
    return ldl_pd([[M[i][j] for j in keep] for i in keep])


cands = []
for wt in [Fr(1, 30), Fr(1, 10), Fr(1), Fr(10), Fr(30)]:
    for (u0, u1) in [(4, 0), (5, 0), (3, 1), (2, 2), (1, 3), (2, 3)]:
        for wq in [Fr(1, 10), Fr(1, 30), Fr(10)]:
            cands.append((wt, (u0, u1), wq))
found = []
for wt, (u0, u1), wq in cands:
    rt = reduced_pd(wt * wt, Fr(1), wq * wq, u0, u1, 1)
    rs = reduced_pd(wt * wt, Fr(1), wq * wq, u0, u1, BS)
    if all(rt.values()) and not all(rs.values()):
        found.append((wt, (u0, u1), wq, [p for p, v in rs.items() if not v]))
print(f'{len(cands)} candidate points scanned with the reduced forms; {len(found)} with a failing scaled part '
      f'[{time.time() - t0:.1f}s]', flush=True)
seen = set()
nconf = 0
for wt, (u0, u1), wq, parts in found:
    key = ((u0, u1), wt in (Fr(1), wq))
    if key in seen:
        continue
    seen.add(key)
    w = [wt] + [Fr(1)] * u0 + [wq] * u1
    pt, ps = block_pd(w, 1), block_pd(w, BS)
    good = pt and not ps
    nconf += good
    ok &= good
    print(f'  w = {[str(x) for x in w]}: failing scaled parts {parts}; j1lib gauge-fixed block PD: true {pt}, scaled '
          f'{ps} -> control {"confirmed" if good else "NOT confirmed"} [{time.time() - t0:.1f}s]', flush=True)
    if nconf >= 6:
        break
ok &= nconf >= 3
tv = [f for f in found if f[0] in (Fr(1), f[2]) and f[1][1] > 0]
print(f'  points with all weights taking two values among the failures: {len(tv)}')
ok &= len(tv) > 0
print('PASS' if ok else 'FAIL')
