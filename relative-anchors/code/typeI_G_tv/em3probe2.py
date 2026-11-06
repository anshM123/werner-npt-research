"""Exploration (no claim rests on it): the cost of one slice computed in the orthant variables and in the variables
of the closed form; both must give the same reduced fraction.  python em3probe2.py s|1 c1 c2 ..."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import em3lib as L

ka = sys.argv[1]
cs = [int(x) for x in sys.argv[2:]]
t0 = time.time()
A, (names, Ao) = L.build_A(ka, orig=True)
OS = L.OrigSlicer(ka, names, Ao)
print('built [%.1fs]' % (time.time() - t0), flush=True)
ok = True
for c in cs:
    t1 = time.time()
    d2 = OS.det(c)
    t2 = time.time()
    d1 = L.slice_det(A, c)
    t3 = time.time()
    same = (d1.n == d2.n) and (d1.d == d2.d)
    ok &= same
    print('x1 = %d: closed-form variables %.1fs, orthant variables %.1fs; equal %s; num %d terms, den %d terms'
          % (c, t2 - t1, t3 - t2, same, len(d2.n), len(d2.d)), flush=True)
print('PASS' if ok else 'FAIL')
