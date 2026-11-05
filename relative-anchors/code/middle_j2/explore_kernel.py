"""Block sizes and kernels of the three families of (R_inf)_2 at k = 5, 6, 7 (exact rational arithmetic).
For each family and weight vector: the size, the kernel vectors p (x) y lie in the kernel, the nullity of the block
equals the number of these vectors (so the kernel is exactly their span), and the block with one coordinate (1@r0, y)
deleted per column y is positive definite (exact LDL).  Weights: family 1 arbitrary on T and two U-levels; family 2
arbitrary on t and two U-levels; family 3 two levels.  Usage: python explore_kernel.py [kmax]"""
import sys, time, random
from fractions import Fraction as Fr
from j1lib import J1
from r2lib import states_T1, states_T2, states_L1, kernel_vectors, rank_exact, ldl_pd

t0 = time.time()
kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
random.seed(5)
ok = True


def run(fam, k, w):
    global ok
    J = J1(w)
    if fam == 1:
        st = states_T2(k, 0, 1)
    elif fam == 2:
        st = states_T1(k, 0)
    else:
        st = states_L1(k)
    M = J.matrix(st)
    cols, kv = kernel_vectors(st, w)
    inker = all(all(sum(M[i][j] * v[j] for j in range(len(st))) == 0 for i in range(len(st))) for v in kv)
    null = len(st) - rank_exact(M)
    # gauge: delete (1@r0, y), r0 = the last copy (a U-copy)
    r0 = k - 1
    keep = [i for i, (a, y) in enumerate(st) if not (a is not None and a[0] == 1 and a[1] == r0)]
    pd = ldl_pd([[M[i][j] for j in keep] for i in keep])
    good = inker and null == len(kv) and pd
    ok &= good
    print(f'family {fam} k={k} w={[str(x) for x in w]}: size {len(st)}, kernel vectors {len(kv)} in kernel: {inker}, '
          f'nullity {null}, gauge-fixed PD: {pd} [{time.time() - t0:.1f}s]', flush=True)


for k in range(5, kmax + 1):
    for trial in range(2):
        a = Fr(random.randint(1, 9), random.randint(1, 9))
        b = Fr(random.randint(1, 9), random.randint(1, 9))
        q1 = Fr(1)
        q2 = Fr(random.randint(1, 9), random.randint(1, 9))
        u1 = random.randint(1, k - 3)
        # family 1: copies 0, 1 = T with weights a, b; U-levels
        w1 = [a, b] + [q1] * u1 + [q2] * (k - 2 - u1)
        run(1, k, w1)
        u1 = random.randint(1, k - 2)
        w2 = [a] + [q1] * u1 + [q2] * (k - 1 - u1)
        run(2, k, w2)
        u1 = random.randint(1, k - 1)
        w3 = [q1] * u1 + [q2] * (k - u1)
        if k <= 6 or trial == 0:
            run(3, k, w3)
print('PASS' if ok else 'FAIL')
