"""A check of Theorem L13 (the letter-1 block of (R_inf)_2 at two weight values, smaller level of 2..13 copies),
independent of its reduction and certificates: the whole letter-1 block of (R_inf)_2 is built from Lemma B by the fast
builder of markcheck.py (equal to j1lib.J1.matrix entry by entry for k <= 6), at random two-value weight vectors with
k copies; the known kernel vectors p (x) |1@X> are annihilated exactly, their number is C(k, 2), and the block without a
set of pivot coordinates of them is positive definite by the exact test pd_cert.
python r2check.py K NPTS SEED"""
import sys, random, time
from fractions import Fraction as Fr
from markcheck import Fast, kernel_ok
from crosscheck import pd_cert, pivots_of, known_kernel
from math import comb

k, npts, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
rng = random.Random(seed)
amps = [Fr(1, 3), Fr(2, 5), Fr(3, 4), Fr(5, 6), Fr(4, 3), Fr(3, 2), Fr(7, 4), Fr(5, 2), Fr(3)]
ok = True
T0 = time.time()
for pt in range(npts):
    b = rng.randint(2, min(13, k // 2))
    lev = [1] * b + [0] * (k - b)
    rng.shuffle(lev)
    a0, a1 = rng.sample(amps, 2)
    w = [a1 if l else a0 for l in lev]
    t0 = time.time()
    states, M = Fast(w).block(2)
    n = len(states)
    vecs = known_kernel(states, k, 2, w)
    kern = kernel_ok(M, vecs)
    cnt = len(vecs) == comb(k, 2)
    piv = set(pivots_of(vecs))
    keep = [i for i in range(n) if i not in piv]
    proved, info = pd_cert([[M[i][jj] for jj in keep] for i in keep])
    good = kern and cnt and proved
    ok &= good
    print(f'k={k} point {pt}: smaller level {b} copies, weights {a0}, {a1} (q = {(a1 / a0) ** 2}): {n} states, '
          f'{len(vecs)} kernel vectors annihilated {kern}, PD complement {proved} ({info}): {good} '
          f'[{time.time() - t0:.1f}s]', flush=True)
print(f'total {time.time() - T0:.0f}s')
print('PASS' if ok else 'FAIL')
