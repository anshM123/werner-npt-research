"""Negative controls: with the B-part (everything that contains an E^{-1}) multiplied by 1 + 1/1000, the block is
indefinite at explicit rational points next to the corner where one copy of T carries almost all the weight, while
the true block is positive definite there.  At the same points the scaled analogue of the pair/row inequalities
(F_il < 0 and a_l F_ii + (x_T - a_i) F_il > 0 for all pairs) fails.  For k <= 5 the scaled block is also taken
directly from ylib as A - (1001/1000) B (units 2^{k-j}) and compared with the scaled closed form."""
import time
import flint
from djlib import block, eliminate_vacuum, states
from ylib import Builder, point_ring
Q = flint.fmpq
t0 = time.time()
lam = Q(1001, 1000)


def pivots(M):
    n = len(M); A = [row[:] for row in M]; out = []
    for c in range(n):
        pv = A[c][c]; out.append(pv)
        if pv <= 0:
            return out
        for i in range(c + 1, n):
            f = A[i][c] / pv
            for l in range(c, n):
                A[i][l] -= f * A[c][l]
    return out


def ineqs_hold(aT, M):
    F = eliminate_vacuum(M)
    n = len(aT); xT = sum(aT, Q(0))
    if not all(F[i][l] < 0 for i in range(n) for l in range(n) if i != l):
        return False
    return all(aT[l] * F[i][i] + (xT - aT[i]) * F[i][l] > 0 for i in range(n) for l in range(n) if i != l)


ok = True
cases = [(1, 0, 10 ** 4), (1, 1, 10 ** 4), (1, 3, 10 ** 4), (2, 0, 10 ** 4), (2, 1, 10 ** 4), (2, 3, 10 ** 5),
         (3, 0, 10 ** 4), (3, 1, 10 ** 4), (3, 3, 10 ** 5), (5, 0, 10 ** 5), (5, 1, 10 ** 5), (5, 3, 10 ** 5)]
for j, N, big in cases:
    s = [Q(big)] + [1 + Q(i, 10) for i in range(j)] + [Q(1, 2) + Q(i, 5) for i in range(N)]
    S = sum(s, Q(0)); x = [v / S for v in s]
    aT, xU = x[:j + 1], x[j + 1:]
    M = block(aT, xU)
    Ms = block(aT, xU, bscale=lam)
    pv = pivots(M); pvs = pivots(Ms)
    true_pd = all(v > 0 for v in pv) and len(pv) == len(M)
    scaled_fails = any(v <= 0 for v in pvs)
    holds_true = ineqs_hold(aT, M)
    holds_scaled = ineqs_hold(aT, Ms)
    good = true_pd and scaled_fails and holds_true and not holds_scaled
    k = j + 1 + N
    extra = ''
    if k <= 5:
        Bd = Builder(point_ring(s))
        sts = states(k, list(range(j + 1)))
        unit = Q(2) ** (k - j)
        Ys = [[unit * (Bd.A(p[0], p[1], q[0], q[1]) - lam * Bd.B(p[0], p[1], q[0], q[1])) for q in sts] for p in sts]
        same = all(Ys[i][l] == Ms[i][l] for i in range(j + 2) for l in range(j + 2))
        good &= same
        extra = f'; ylib A - (1001/1000) B equals the scaled closed form: {same}'
    ok &= good
    print(f'j = {j}, k = {k}, s = ({big}, ...): true block positive definite {true_pd}; scaled block has pivot '
          f'{float(pvs[-1]):.3e} <= 0 at position {len(pvs)}: {scaled_fails}; pair/row inequalities true {holds_true}, '
          f'scaled {holds_scaled}{extra} [{time.time() - t0:.1f}s]', flush=True)
print('PASS' if ok else 'FAIL')
