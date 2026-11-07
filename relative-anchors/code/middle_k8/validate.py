"""FastJ1.matrix against j1lib.J1.matrix, entry by entry, on every torus group of every (R_inf)_j.
python validate.py K NPTS SEED"""
import sys, random, time
from fractions import Fraction as Fr
from j1lib import J1
from fastj1 import FastJ1
from crosscheck import rinf_states, label

k, npts, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
rng = random.Random(seed)
amps = [Fr(1, 3), Fr(2, 5), Fr(3, 4), Fr(5, 6), Fr(4, 3), Fr(3, 2), Fr(7, 4), Fr(5, 2), Fr(3)]
ok = True
t0 = time.time()
ngroups = nent = 0
for pt in range(npts):
    w = [rng.choice(amps) for _ in range(k)]
    for bs in (Fr(1), Fr(1001, 1000)):
        J, F = J1(w, bscale=bs), FastJ1(w, bscale=bs)
        for j in range(-1, k + 1):
            groups = {}
            for s in rinf_states(k, j):
                groups.setdefault(label(s, k), []).append(s)
            for lab, gs in groups.items():
                A = J.matrix(gs)
                B = F.matrix(gs)
                n = len(gs)
                same = all(A[p][q] == B[p][q] for p in range(n) for q in range(n))
                ngroups += 1
                nent += n * n
                if not same:
                    ok = False
                    print('DIFFERS', k, w, bs, j, lab, flush=True)
    print(f'point {pt}: weights {[str(x) for x in w]} done [{time.time() - t0:.1f}s]', flush=True)
print(f'{ngroups} groups, {nent} entries compared')
print('PASS' if ok else 'FAIL')
