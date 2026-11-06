"""A check of Theorems N5, N6, E and G3 (letter-1 blocks with one marked copy, parts of irrep index i >= 5) at members
beyond every fixed-k certificate (56 <= k <= KMAX), independent of the family certificates: the numeric reduced form Q_i
of mlib.py (checked against projections of the Lemma B block for k <= 9 by projcheck.py) at the member (k, j, i) and a
rational sigma is positive semidefinite with the expected nullity: exact rank over Q, and a principal submatrix of size
n - nullity that is positive definite by exact symmetric elimination.  Expected nullity (the nonzero kappa_y0):
0 if q0 = -1 or t = -2, 1 if q0 = 0 or t = -1, else 2.
Regions sampled: 'n56' = [q0 >= 4, t >= 4, i in {5, 6}]; 'eq' = [q0 <= 3, t >= 4, i >= 6]; 'et' = [q0 >= 4, t <= 3,
i >= 5]; 'g7' = [q0 >= 4, t >= 4, i >= 7] (the whole form, whose two-row part is not proved; a check only).
python mfcheck.py REGION NPTS SEED KMAX"""
import sys, os, random, time
from fractions import Fraction as Fr
import flint
import mlib

region, npts, seed, kmax = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
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


def sample():
    if region == 'n56':
        return rng.randint(4, 60), rng.randint(4, 60), rng.choice([5, 6])
    if region == 'eq':
        return rng.randint(-1, 3), rng.randint(4, 60), rng.randint(6, 40)
    if region == 'et':
        return rng.randint(4, 60), rng.randint(-2, 3), rng.randint(5, 40)
    return rng.randint(4, 40), rng.randint(4, 40), rng.randint(7, 30)


done = 0
tries = 0
while done < npts and tries < 10 ** 6:
    tries += 1
    q0, t, i = sample()
    k = q0 + t + 2 * i + 1
    j = q0 + i
    if not (56 <= k <= kmax and 9 <= j <= k - 9):
        continue
    null = 0 if (q0 == -1 or t == -2) else (1 if (q0 == 0 or t == -1) else 2)
    sg = rng.choice(sigmas)
    names, norms, Q = mlib.reduced(mlib.Calc(k, j, i, sigma=sg))
    n = len(names)
    A = flint.fmpq_mat(n, n, [Q[a][c] for a in range(n) for c in range(n)])
    rank = A.rank()
    good = False
    if rank == n - null:
        # a principal submatrix of size n - null with nonzero determinant, tested for positive definiteness
        import itertools
        for drop in itertools.combinations(range(n), null):
            keep = [a for a in range(n) if a not in drop]
            sub = flint.fmpq_mat(len(keep), len(keep), [Q[a][c] for a in keep for c in keep])
            if len(keep) == 0 or sub.det() != 0:
                good = pd([[Q[a][c] for c in keep] for a in keep])
                break
    res = rank == n - null and good
    ok &= res
    done += 1
    print(f'{region} k={k} j={j} i={i} (q0={q0}, t={t}) sigma={sg}: {n} twisted vectors, rank {rank} (nullity expected '
          f'{null}), PD complement {good}: {res} [{time.time() - T0:.1f}s]', flush=True)
print('PASS' if ok else 'FAIL')
