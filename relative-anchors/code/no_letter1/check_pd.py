"""Exact checks of the conclusions at rational weight vectors (closed form of djlib, exact arithmetic):
  the block M is positive definite (all pivots of a symmetric elimination positive), F (vacuum eliminated) has
  negative off-diagonal entries, every row satisfies F_ii > sum_{l != i} |F_il|, and every pair satisfies
  a_l F_ii + (x_T - a_i) F_il > 0.  Weight vectors: random, near the corner where one copy of T carries almost all
  the weight, with tiny or huge U-weights, and with many copies in U; |T| = j + 1 for j = 1, ..., 7, |U| up to 12.
Also, for k <= 6, positive definiteness of the block taken directly from ylib at a few of these weight vectors."""
import random, time
import flint
from djlib import block, eliminate_vacuum, states
from ylib import Builder, point_ring
Q = flint.fmpq
random.seed(1105)
t0 = time.time()


def pivots(M):
    n = len(M); A = [row[:] for row in M]; out = []
    for c in range(n):
        pv = A[c][c]; out.append(pv)
        if pv == 0:
            return out
        for i in range(c + 1, n):
            f = A[i][c] / pv
            if f != 0:
                for l in range(c, n):
                    A[i][l] -= f * A[c][l]
    return out


def weights(n, N, kind):
    if kind == 'random':
        w = [Q(random.randint(1, 50), random.randint(1, 50)) for _ in range(n + N)]
    elif kind == 'corner':
        w = [Q(random.randint(1, 9), random.randint(1, 9)) for _ in range(n + N)]
        w[0] = w[0] * random.choice([10 ** 3, 10 ** 5, 10 ** 7])
    elif kind == 'tinyU':
        w = [Q(random.randint(1, 9)) for _ in range(n)] + [Q(1, random.randint(10 ** 3, 10 ** 6)) for _ in range(N)]
    elif kind == 'hugeU':
        w = [Q(1, random.randint(1, 9)) for _ in range(n)] + [Q(random.randint(10 ** 2, 10 ** 4)) for _ in range(N)]
    else:  # spread
        w = [Q(random.randint(1, 9), 10 ** random.randint(0, 6)) for _ in range(n + N)]
    S = sum(w, Q(0))
    return [v / S for v in w]


ok = True; cnt = 0
for j in range(1, 8):
    n = j + 1
    for N in (0, 1, 2, 5, 12):
        for kind in ('random', 'corner', 'tinyU', 'hugeU', 'spread'):
            if N == 0 and kind in ('tinyU', 'hugeU'):
                continue
            w = weights(n, N, kind)
            aT, xU = w[:n], w[n:]
            M = block(aT, xU)
            pd = all(pv > 0 for pv in pivots(M))
            F = eliminate_vacuum(M)
            xT = sum(aT, Q(0))
            offneg = all(F[i][l] < 0 for i in range(n) for l in range(n) if i != l)
            rows = all(F[i][i] > sum((-F[i][l] for l in range(n) if l != i), Q(0)) for i in range(n))
            pairs = all(aT[l] * F[i][i] + (xT - aT[i]) * F[i][l] > 0 for i in range(n) for l in range(n) if i != l)
            good = pd and offneg and rows and pairs
            ok &= good; cnt += 1
            if not good:
                print('FAILED at', j, N, kind, pd, offneg, rows, pairs, flush=True)
    print(f'j = {j}: all checks passed so far: {ok} [{time.time() - t0:.1f}s]', flush=True)
# directly from ylib, small k
for (k, j) in [(3, 1), (4, 1), (4, 2), (5, 2), (5, 3), (6, 3), (6, 4), (6, 5)]:
    for kind in ('random', 'corner', 'spread'):
        w = weights(j + 1, k - j - 1, kind)
        s = w  # squared weights, |s| = 1
        Bd = Builder(point_ring(s))
        T = list(range(j + 1))
        sts = states(k, T)
        Y = [[Bd.Y(p, q) for q in sts] for p in sts]
        pd = all(pv > 0 for pv in pivots(Y))
        ok &= pd; cnt += 1
        print(f'ylib k = {k}, j = {j}, {kind}: block positive definite {pd} [{time.time() - t0:.1f}s]', flush=True)
print(cnt, 'weight vectors'); print('PASS' if ok else 'FAIL')
