"""Check: the letter-2 orbit calculus (orbit2.py, Num2) reproduces eqlib2 (explicit vectors; itself checked against
the independent builder for k = 3, 4 and used by run2 for k <= 14) exactly, for every case I form with k <= KMAX."""
import sys, time
from fractions import Fraction
from orbit2 import Num2, reduced2
from eqlib2 import caseI_vectors, form
KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
ok = True; n = 0; t0 = time.time()
for k in range(2, KMAX + 1):
    for j in range(0, k + 1):
        for tau in range(1, min(k, j + 1) + 1):
            for i in range(0, (k - tau) // 2 + 1):
                for mode in (('sym', 'anti') if tau >= 2 else ('sym',)):
                    vecs = caseI_vectors(k, j, tau, i, mode)
                    names2, G2, Q2 = reduced2(Num2(k, j, tau, i), mode)
                    names1 = [v[0] for v in vecs]
                    if not vecs and not names2:
                        continue
                    G1, Q1 = form(k, j, vecs) if vecs else ([], [])
                    conv = lambda M: [[Fraction(int(x.p), int(x.q)) for x in row] for row in M]
                    n += 1
                    if names1 != names2 or conv(G1) != G2 or conv(Q1) != Q2:
                        ok = False
                        print('MISMATCH', k, j, tau, i, mode, names1, names2)
    print(f'k={k} done ({time.time()-t0:.0f}s)', flush=True)
print(n, 'forms compared;', 'PASS' if ok else 'FAIL')
