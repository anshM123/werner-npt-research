"""Negative controls inside the regimes of the extension (b + x, b), 14 <= b <= 40 (part (b) of neg.py at larger
sizes): at (u0, u1) = (30, 14), (25, 20), (41, 27), (50, 33), (41, 40), every W-block of the reduced form is positive
definite for the true B and some W-block is not for B -> (1 + 1/1000)B, at an explicit q.  PASS iff all of this
holds."""
import time
from fractions import Fraction as Fr
import flint
from fastj1 import FastJ1
from j1lib import mu
from r2lib import states_L1
import r2v

t0 = time.time()
B1 = Fr(1001, 1000)
ok = True


def fq(x):
    return flint.fmpq(x.numerator, x.denominator)


def tofr(x):
    return Fr(int(x.p), int(x.q))


def wblocks(u0, u1, q, b):
    """dict family -> W-block of the reduced S_b (r2v.reduce_scaled)"""
    out = {}
    ctx = r2v.Ctx(Fr(1), q, Fr(u0), Fr(u1))
    for F in r2v.families():
        bases = r2v.choose(F, (u0, u1))
        if not bases[1]:
            continue
        S = r2v.reduce_scaled(F, ctx, bases, b)
        nK = len(bases[0])
        out[F.name] = (S, [row[nK:] for row in S[nK:]])
    return out


# (b) inside the regimes
for (u0, u1) in ((30, 14), (25, 20), (41, 27), (50, 33), (41, 40)):
    found = None
    for q in (Fr(1, 10 ** 4), Fr(10 ** 4), Fr(1, 100), Fr(100), Fr(1, 10 ** 6), Fr(10 ** 6)):
        t = wblocks(u0, u1, q, Fr(1))
        if not all(r2v.pd_exact(v[1]) for v in t.values()):
            ok = False
            print(f'(b) ({u0}, {u1}) q = {q}: a W-block is NOT PD for the true B (contradiction)')
            break
        s = wblocks(u0, u1, q, B1)
        bad = [nm for nm, v in s.items() if not r2v.pd_exact(v[1])]
        if bad:
            found = (q, bad)
            break
    ok &= found is not None
    print(f'(b) ({u0}, {u1}): all W-blocks PD at b = 1; at b = 1001/1000, q = {found[0] if found else None}: '
          f'not PD for {found[1] if found else None} [{time.time() - t0:.1f}s]', flush=True)
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
