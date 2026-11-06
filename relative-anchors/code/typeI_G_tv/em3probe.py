"""Exploration for the determinant route (no claim rests on it): entry sizes of A, the factors of the entry
denominators, the exponents of D*, the degree bounds, and the cost and the denominator exponents of a few slices.
python em3probe.py s|1 c1 c2 ..."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import em3lib as L

ka = sys.argv[1]
cs = [int(x) for x in sys.argv[2:]]
t0 = time.time()
A = L.build_A(ka)
factors, mult, estar = L.factor_data(A)
print('factors of the entry denominators (%d) [%.1fs]:' % (len(factors), time.time() - t0))
for k, f in enumerate(factors):
    tab = [[mult.get((k, a, b), 0) for b in range(3)] for a in range(3)]
    print('  f%d: %d terms, degrees %s, positive coefficients and a monomial free of x1, x2, r: %s; '
          'multiplicities %s; e_f = %d' % (k, len(f), f.degrees(), L.one_signed_pos(f, ('x1', 'x2', 'r')), tab,
                                          estar[k]))
for iv, v in enumerate(L.VN):
    B, dD = L.degree_bound(A, factors, estar, iv)
    print('  %s: deg D* = %d, bound B = %d' % (v, dD, B))
print('num degrees of the entries:', [A[i][j].n.degrees() for i in range(3) for j in range(i, 3)])
print('den degrees of the entries:', [A[i][j].d.degrees() for i in range(3) for j in range(i, 3)])
for c in cs:
    t1 = time.time()
    dt = L.slice_det(A, c)
    t2 = time.time()
    obs = []
    for f in factors:
        fs = L.spec(f, c)
        if fs.is_constant():
            obs.append('-')
            continue
        obs.append(L.multiplicity(fs, dt.d))
    print('slice x1 = %d: det num %d terms, degrees %s; den %d terms, degrees %s; det %.1fs; observed exponents %s; '
          'num coefficients all positive: %s [%.1fs]'
          % (c, len(dt.n), dt.n.degrees(), len(dt.d), dt.d.degrees(), t2 - t1, obs,
             all(x > 0 for x in dt.n.coeffs()), time.time() - t1), flush=True)
print('done [%.1fs]' % (time.time() - t0))
print('PASS')
