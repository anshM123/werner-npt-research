"""fastpd.ldl_pd_fast against the sequential exact elimination (the original k5lib.ldl_pd for size <= 24 and the same
loop written out here for any size), on random symmetric rational matrices, positive definite or not: the same
(ok, index of the first nonpositive pivot).  Also the flint matrix-vector product of k5lib.matvec against the sum.
python validate_pd.py NTRIALS SEED"""
import sys, random
from fractions import Fraction as Fr
import k5lib as K
from fastpd import ldl_pd_fast


def ldl_seq(M):
    A = [[Fr(x) for x in row] for row in M]
    n = len(A)
    for c in range(n):
        p = A[c][c]
        if p <= 0:
            return False, c
        for i in range(c + 1, n):
            if A[i][c] != 0:
                f = A[i][c] / p
                for l in range(c + 1, n):
                    if A[c][l] != 0:
                        A[i][l] -= f * A[c][l]
    return True, None


ntr, seed = int(sys.argv[1]), int(sys.argv[2])
rng = random.Random(seed)
ok = True
npd = 0
for trial in range(ntr):
    n = rng.choice([5, 23, 25, 40, 61, 90, 130])
    X = [[Fr(rng.randint(-5, 5), rng.randint(1, 4)) for _ in range(n)] for _ in range(n)]
    M = [[sum(X[i][t] * X[j][t] for t in range(n)) for j in range(n)] for i in range(n)]
    if trial % 3 == 1:
        sh = Fr(rng.randint(1, 40))
        for i in range(n):
            M[i][i] -= sh
    a, b = ldl_seq(M), ldl_pd_fast(M)
    npd += a[0]
    if a != b:
        ok = False
        print('DIFFERS', n, a, b)
    v = [Fr(rng.randint(-3, 3), rng.randint(1, 3)) for _ in range(n)]
    if K.matvec(M, v) != [sum(M[i][l] * v[l] for l in range(n)) for i in range(n)]:
        ok = False
        print('MATVEC DIFFERS', n)
print(f'{ntr} matrices ({npd} positive definite): same result and same first nonpositive pivot')
print('PASS' if ok else 'FAIL')
