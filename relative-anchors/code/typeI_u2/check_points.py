"""Exact point checks (sanity; Theorem C covers all cases): the principal block of the type I block (cf.TypeI) on
{P_{r,h}: r, h in U} u {Q} has all LDL^T pivots positive (Fractions), for random rational weights with u = 2..5,
m = 1..3, including widely spread weights (ratios up to 10^4).  Prints PASS."""
import sys, os, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from cf import TypeI

random.seed(13)
ok = True
t0 = time.time()
cnt = 0
for u in (2, 3, 4, 5):
    for m in (1, 2, 3):
        if u + m > 7:
            continue
        for trial in range(3):
            w = [Fr(random.randint(1, 40), random.randint(1, 40)) for _ in range(u + m)]
            if trial == 1:
                w[0] *= 100          # one T-weight large
            if trial == 2:
                w[m] *= 100          # one U-weight large
            T = list(range(m))
            ti = TypeI(w, T)
            M = ti.build()
            keep = [ti.idx[('P', r, h)] for r in ti.U for h in ti.U] + [ti.idx[('Q',)]]
            A = [[M[i][j] for j in keep] for i in keep]
            n = len(A)
            piv_ok = True
            for i in range(n):
                if A[i][i] <= 0:
                    piv_ok = False
                    break
                for j in range(i + 1, n):
                    f = A[j][i] / A[i][i]
                    if f != 0:
                        for k in range(i, n):
                            A[j][k] -= f * A[i][k]
            ok &= piv_ok
            cnt += 1
            print(f"u={u} m={m} trial {trial}: all {n} pivots positive: {piv_ok}", flush=True)
print(f"{cnt} blocks, time {time.time() - t0:.1f}s")
print("PASS" if ok else "FAIL")
