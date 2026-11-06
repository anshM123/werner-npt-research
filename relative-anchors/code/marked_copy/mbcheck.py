"""A check of Theorems N1-N4 at members beyond every fixed-k certificate (41 <= k <= 60), independent of the family
certificates: the numeric reduced form Q_i of mlib.py (checked against projections of the Lemma B block for k <= 9 by
projcheck.py) at the member (k, j, i) and a rational sigma is positive semidefinite with nullity exactly 2
(= dim of K_j on V_i): exact rank over Q, and a principal submatrix of size n - 2 that is positive definite by exact
symmetric elimination.
python mbcheck.py NPTS SEED"""
import sys, os, random, time
from fractions import Fraction as Fr
import flint
import mlib

npts, seed = int(sys.argv[1]), int(sys.argv[2])
rng = random.Random(seed)
sigmas = [Fr(1, 10 ** 6), Fr(1, 97), Fr(3, 7), Fr(1), Fr(5, 2), Fr(64), Fr(10 ** 6)]
ok = True
T0 = time.time()


def pd(M):
    n = len(M)
    A = [row[:] for row in M]
    for c in range(n):
        if A[c][c] <= 0:
            return False
        for i in range(c + 1, n):
            if A[i][c] != 0:
                f = A[i][c] / A[c][c]
                for j in range(c + 1, n):
                    A[i][j] -= f * A[c][j]
    return True


done = 0
while done < npts:
    b = rng.choice([1, 2, 3, 4])
    q0 = rng.randint(4, 40)
    t = rng.randint(4, 40)
    k = q0 + t + 2 * b + 1
    if not 41 <= k <= 60:
        continue
    j = q0 + b
    if not 9 <= j <= k - 9:
        continue
    sg = rng.choice(sigmas)
    names, norms, Q = mlib.reduced(mlib.Calc(k, j, b, sigma=sg))
    n = len(names)
    A = flint.fmpq_mat(n, n, [Q[a][c] for a in range(n) for c in range(n)])
    rank = A.rank()
    found = False
    for d1 in range(n):
        for d2 in range(d1 + 1, n):
            keep = [a for a in range(n) if a not in (d1, d2)]
            sub = flint.fmpq_mat(n - 2, n - 2, [Q[a][c] for a in keep for c in keep])
            if sub.det() != 0:
                found = True
                good = pd([[Q[a][c] for c in keep] for a in keep])
                break
        if found:
            break
    res = found and good and rank == n - 2
    ok &= res
    done += 1
    print(f'k={k} j={j} i={b} (q0={q0}, t={t}) sigma={sg}: {n} twisted vectors, rank {rank}, PD complement {found and good}: '
          f'{res} [{time.time() - T0:.1f}s]', flush=True)
print('PASS' if ok else 'FAIL')
