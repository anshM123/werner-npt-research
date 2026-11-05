"""Exact checks of Conjecture L1 (NOTES.md, Section 4) for the letter-1 block of (R_inf)_1, k >= 4.
States P_{r,c} = (1@r, 1@c), Q_{cc'} = (vac, 1@{c,c'}); rescaled coordinates of ylib (P_{r,c} divided by w_r w_c,
Q_{cc'} by w_c w_c').  Stage 0: eliminate the vacuum states (their block is positive definite by Lemma G).  Then for
the column groups c in a given order: the diagonal block of the current Schur complement on the group
{P_{r,c} : r} must be a Laplacian (zero row sums, every off-diagonal entry < 0); then the group is eliminated in the
gauge (its coordinate with the smallest r is dropped; the rest of a connected Laplacian is positive definite).  If
every stage passes, the block is >= 0 with kernel exactly p (x) (letter-1 part of H_1) at that weight vector.
Exact arithmetic (Fractions), rational weights.  Also a negative control: with B scaled by 1 + 1/1000 the gauge-fixed
block (one coordinate (1@r0, 1@c) dropped per column) is not positive definite at an explicit rational point.
Usage: python check_l1.py [kmax] [points_per_k]"""
import sys, random
from fractions import Fraction as Fr
from j1lib import J1, block_states, mu

kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
npts = int(sys.argv[2]) if len(sys.argv) > 2 else 4
random.seed(11)


def rescaled(J, st):
    M = J.matrix(st)
    m = [mu(z, J.w) for z in st]
    return [[m[i] * M[i][j] * m[j] for j in range(len(st))] for i in range(len(st))]


def schur(M, keep, elim):
    n = len(elim)
    A = [[M[i][j] for j in elim] + [M[i][j] for j in keep] for i in elim]
    for c in range(n):
        p = A[c][c]
        if p <= 0:
            return None
        A[c] = [v / p for v in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[c])]
    X = [row[n:] for row in A]
    return [[M[keep[a]][keep[b]] - sum(M[keep[a]][elim[l]] * X[l][b] for l in range(n) if X[l][b] != 0)
             for b in range(len(keep))] for a in range(len(keep))]


def run(w, order, bscale=1):
    k = len(w)
    J = J1(w, bscale=bscale)
    st = block_states(k, ())
    M = rescaled(J, st)
    P = [i for i, (a, y) in enumerate(st) if a is not None]
    Q = [i for i, (a, y) in enumerate(st) if a is None]
    S = schur(M, P, Q)
    if S is None:
        return False, 'vacuum block not PD'
    cur = P[:]
    col = {i: st[i][1].index(1) for i in P}
    for c in order:
        g = [a for a, i in enumerate(cur) if col[i] == c]
        D = [[S[a][b] for b in g] for a in g]
        if any(sum(row) != 0 for row in D):
            return False, f'row sums of group {c} not zero'
        if any(D[a][b] >= 0 for a in range(len(g)) for b in range(len(g)) if a != b):
            return False, f'group {c} not a Laplacian'
        el = g[1:]
        rest = [a for a in range(len(cur)) if a not in g]
        if rest:
            S2 = schur(S, rest, el)
            if S2 is None:
                return False, 'pivot'
            S = S2
        cur = [cur[a] for a in rest]
    return True, 'ok'


ok = True
cnt = 0
for k in range(4, kmax + 1):
    for i in range(npts):
        mode = i % 4
        if mode == 0:
            w = [Fr(random.randint(1, 30), random.randint(1, 30)) for _ in range(k)]
        elif mode == 1:
            w = [Fr(2) ** random.randint(-5, 5) * Fr(random.randint(1, 3)) for _ in range(k)]
        elif mode == 2:
            w = [Fr(random.randint(1, 9)) for _ in range(k)]; w[random.randrange(k)] = Fr(1, 200)
        else:
            w = [Fr(1, random.randint(1, 9)) for _ in range(k)]; w[random.randrange(k)] = Fr(40)
        for order in (list(range(k)), random.sample(range(k), k)):
            good, msg = run(w, order)
            cnt += 1
            if not good:
                print(f'FAIL k={k} w={[str(v) for v in w]} order={order}: {msg}', flush=True)
            ok &= good
    print(f'k = {k}: done, all stages pass so far: {ok}', flush=True)
print(f'{cnt} (weight vector, order) pairs: every stage is a Laplacian: {ok}')
# negative control: with B -> (1 + 1/1000) B the gauge-fixed block (one coordinate (1@r0, 1@c) dropped per column)
# is not positive definite at an explicit rational point (exact LDL), while the true one is.


def ldl_pd(A):
    A = [row[:] for row in A]
    n = len(A)
    for c in range(n):
        if A[c][c] <= 0:
            return False
        for i in range(c + 1, n):
            if A[i][c] != 0:
                f = A[i][c] / A[c][c]
                for j in range(c + 1, n):
                    A[i][j] -= f * A[c][j]
    return True


found = None
fams = [lambda a: [Fr(1), Fr(1), Fr(1), Fr(1, a)], lambda a: [Fr(1), Fr(1, a), Fr(1, a + 1), Fr(1, a + 2)],
        lambda a: [Fr(a), Fr(1), Fr(2), Fr(3)], lambda a: [Fr(1), Fr(1), Fr(1, a), Fr(1, a)]]
cands = [f(a) for a in (2, 3, 5, 8, 13, 20, 30, 50, 80, 130, 200, 300) for f in fams]
for w in cands:
    st = block_states(4, ())
    keep = [i for i, (aa, y) in enumerate(st) if not (aa is not None and aa[1] == 0)]
    Ms = rescaled(J1(w, bscale=Fr(1001, 1000)), st)
    M0 = rescaled(J1(w), st)
    if not ldl_pd([[Ms[i][j] for j in keep] for i in keep]) and ldl_pd([[M0[i][j] for j in keep] for i in keep]):
        found = [str(v) for v in w]
        break
print(f'negative control: B -> (1 + 1/1000) B makes the gauge-fixed block not PD at k = 4, w = {found}')
ok &= found is not None
print('PASS' if ok else 'FAIL')
