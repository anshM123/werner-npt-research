"""A check of Corollaries TV and M1 ((R_inf)_{k-1} at two weight values), independent of their closed forms and
certificates: at amplitudes (A, ..., A, 1, ..., 1) (N copies of amplitude A, a rational, the other k - N of amplitude 1),
every torus group of (R_inf)_j for the given j is built from Lemma B (the fast builder FastJ1), its known kernel
vectors (crosscheck.known_kernel) are annihilated exactly, the group without a set of pivot coordinates of them is
positive definite (exact test crosscheck.check_group), and the kernel dimensions add up to |H_j| (+ 1 at j = k - 2,
+ k - 1 at j = k - 1).  As markall.py, with N copies of the other weight instead of one.
python markall2.py K N A J1,J2,..."""
import sys, time
from fractions import Fraction as Fr
from crosscheck import rinf_states, label, known_kernel, check_group, strings
from fastj1 import FastJ1

k, nodd, amp = int(sys.argv[1]), int(sys.argv[2]), Fr(sys.argv[3])
js = [int(x) for x in sys.argv[4].split(',')]
w = [amp] * nodd + [Fr(1)] * (k - nodd)
B = FastJ1(w)
allok = True
T0 = time.time()
for j in js:
    t0 = time.time()
    st = rinf_states(k, j)
    groups = {}
    for s in st:
        groups.setdefault(label(s, k), []).append(s)
    ok = True
    nker = 0
    kinds = {}
    for lab_, gs in groups.items():
        M = B.matrix(gs)
        vecs = known_kernel(gs, k, j, w)
        kern, proved, info = check_group(M, vecs)
        nker += len(vecs)
        if max(lab_) >= 2:
            kind = 'type II'
        elif sum(lab_) == 0:
            kind = 'letter-1'
        elif sum(lab_) <= j:
            kind = 'type I with a letter 1'
        else:
            kind = 'type I without a letter 1'
        n_ok, n_all = kinds.get(kind, (0, 0))
        kinds[kind] = (n_ok + (1 if (kern and proved) else 0), n_all + 1)
        if not (kern and proved):
            ok = False
            print(f'   FAILED group {lab_} (size {len(gs)}): kernel {kern}, PD {proved} ({info})', flush=True)
    expected = len(strings(k, j)) + (1 if j == k - 2 else 0) + (k - 1 if j == k - 1 else 0)
    ok &= nker == expected
    allok &= ok
    desc = '; '.join(f'{kd}: {a}/{b}' for kd, (a, b) in sorted(kinds.items()))
    print(f'k={k}, {nodd} copies of amplitude {amp}, (R_inf)_{j}: {len(st)} states in {len(groups)} torus groups '
          f'({desc} passed); kernel dimension {nker} (expected {expected}): {ok} [{time.time() - t0:.1f}s]',
          flush=True)
print(f'total {time.time() - T0:.0f}s')
print('PASS' if allok else 'FAIL')
