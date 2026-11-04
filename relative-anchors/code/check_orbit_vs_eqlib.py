"""Check: the orbit calculus (orbit.py, no enumeration) reproduces eqlib.reduced (explicit vectors, itself
checked against the independent builder for k <= 5) exactly, for k = 2..9, every j = 0..k, every i."""
import time
from fractions import Fraction
from orbit import Num, reduced_orbit
from eqlib import reduced
ok = True
t0 = time.time()
for k in range(2, 10):
    for j in range(0, k + 1):
        for i in range(0, k // 2 + 1):
            n1, G1, Q1, _ = reduced(k, j, i)
            n2, G2, Q2 = reduced_orbit(Num(k, j, i))
            conv = lambda M: [[Fraction(int(x.p), int(x.q)) for x in row] for row in M]
            if n1 != n2 or conv(G1) != G2 or conv(Q1) != Q2:
                ok = False
                print('MISMATCH', k, j, i, n1, n2)
print('PASS' if ok else 'FAIL', f'{time.time()-t0:.1f}s')
