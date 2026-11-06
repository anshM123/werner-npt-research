"""A check of Theorem I1 and Corollary I2 (one marked copy: all of (R_inf)), independent of their reduction and
certificates: at squared weights (a^2, 1, ..., 1) on k copies, every torus group of (R_inf)_j for the given j is built
from Lemma B (the fast builder FastJ1, equal to j1lib.J1 entry by entry for k <= 6); its known kernel vectors
(crosscheck.known_kernel) are annihilated exactly, the group without a set of pivot coordinates of them is positive
definite (exact test crosscheck.pd_cert), and the kernel dimensions add up to |H_j| (+ 1 at j = k - 2, + k - 1 at
j = k - 1).  The groups are classified as letter-1 / type I with a letter 1 / without a letter 1 / type II.
python markall.py K A J1,J2,...   (A = the amplitude of copy 0, a rational)"""
import sys, time
from fractions import Fraction as Fr
from crosscheck import rinf_states, label, known_kernel, check_group, strings
from fastj1 import FastJ1

k, amp, js = int(sys.argv[1]), Fr(sys.argv[2]), [int(x) for x in sys.argv[3].split(',')]
w = [amp] + [Fr(1)] * (k - 1)
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
    print(f'k={k}, amplitude of copy 0 = {amp}, (R_inf)_{j}: {len(st)} states in {len(groups)} torus groups ({desc} '
          f'passed); kernel dimension {nker} (expected {expected}): {ok} [{time.time() - t0:.1f}s]', flush=True)
print(f'total {time.time() - T0:.0f}s')
print('PASS' if allok else 'FAIL')
