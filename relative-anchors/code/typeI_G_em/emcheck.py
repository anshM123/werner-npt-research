"""A check of Corollary EM (type I blocks of (R_inf)_(k-1) whose copies without a letter 2 have equal weight),
independent of its closed forms and certificates: at random weight vectors that are equal on a set U0 of u >= 3 copies
and arbitrary on the other m = k - u >= 2 copies, every torus group of (R_inf)_(k-1) is built from Lemma B (the
reference builder j1lib.J1 for k <= 6, the fast builder FastJ1 for k >= 7, equal to J1 entry by entry for k <= 6);
its known kernel vectors are annihilated exactly, the group without a set of pivot coordinates of them is positive
definite (exact test pd_cert), and the kernel dimensions add up to |H_(k-1)| + (k - 1).  The type I group of
T = complement of U0 is reported separately.
python emcheck.py K NPTS SEED"""
import sys, random, time
from fractions import Fraction as Fr
from crosscheck import rinf_states, label, known_kernel, check_group, strings

k, npts, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
if k <= 6:
    from j1lib import J1 as Builder
else:
    from fastj1 import FastJ1 as Builder
rng = random.Random(seed)
amps = [Fr(1, 10), Fr(1, 3), Fr(2, 5), Fr(3, 4), Fr(5, 6), Fr(1), Fr(4, 3), Fr(3, 2), Fr(7, 4), Fr(5, 2), Fr(3), Fr(7)]
j = k - 1
allok = True
T0 = time.time()
for pt in range(npts):
    u = rng.randint(3, k - 2)
    U0 = sorted(rng.sample(range(k), u))
    om = rng.choice(amps)
    w = [om if r in U0 else rng.choice(amps) for r in range(k)]
    B = Builder(w)
    t0 = time.time()
    st = rinf_states(k, j)
    groups = {}
    for s in st:
        groups.setdefault(label(s, k), []).append(s)
    ok = True
    nker = 0
    typeI = None
    for lab_, gs in groups.items():
        M = B.matrix(gs)
        vecs = known_kernel(gs, k, j, w)
        kern, proved, info = check_group(M, vecs)
        nker += len(vecs)
        if tuple(lab_) == tuple(0 if r in U0 else 1 for r in range(k)):
            typeI = (len(gs), len(vecs), kern and proved)
        if not (kern and proved):
            ok = False
            print(f'   FAILED group {lab_} (size {len(gs)}): kernel {kern}, PD {proved} ({info})', flush=True)
    expected = len(strings(k, j)) + (k - 1)
    ok &= nker == expected and typeI is not None and typeI[2]
    allok &= ok
    print(f'k={k} point {pt}: U = {U0} (u = {u}, m = {k - u}), weights {[str(x) for x in w]}: {len(groups)} torus groups '
          f'of (R_inf)_{j}, each >= 0 with kernel exactly its known vectors; kernel dimension {nker} (expected {expected}); '
          f'the type I group of T = complement of U: {typeI[0]} states, {typeI[1]} kernel vectors, OK {typeI[2]}: {ok} '
          f'[{time.time() - t0:.1f}s]', flush=True)
print(f'total {time.time() - T0:.0f}s')
print('PASS' if allok else 'FAIL')
