"""50-digit exploration of the principal sub-block C on the R- and Q-states (NOTES.md, Section 5); no claim rests on it.

C > 0 if and only if A1 N A1 > (1/2) H A2 H^T (congruence by diag(A1, A2), A1 = I + N1/2, A2 = I + N'/2, then the
Q-states eliminated with (2 A2)^{-1}).  Printed: the largest generalised eigenvalue r of (1/2) H A2 H^T against
A1 N A1 (C > 0 iff r < 1), in the limit sigma -> 0 with two copies of U, and over 300 rational weight vectors (k = 3..7).
Closed forms of tlib.py in mpmath (50 digits); double precision is not reliable here (it gives r > 1 at sigma = 1e-6,
k = 3, a cancellation artefact).
Usage: python explore_rq.py"""
import random
from fractions import Fraction as Fr
import mpmath as mp
from tlib import closed, mul, T

mp.mp.dps = 50


def ratio(xs):
    w = [mp.sqrt(mp.mpf(v.numerator) / v.denominator) for v in xs]
    Cf = closed(w, 0, one=mp.mpf(1))
    A1, A2, H, N = Cf['A1'], Cf['A2'], Cf['H'], Cf['N']
    k = len(xs)
    Nd = [[N[i] if i == j else 0 for j in range(k)] for i in range(k)]
    L = mp.matrix(mul(mul(A1, Nd), A1))
    R = mp.matrix(mul(mul(H, A2), T(H))) * mp.mpf('0.5')
    Li = mp.cholesky(L) ** -1
    return max(mp.eigsy(Li * R * Li.T, eigvals_only=True))


for s in (Fr(1, 10 ** 2), Fr(1, 10 ** 4), Fr(1, 10 ** 6)):
    print(f'k = 3, x = (sigma, 2/5, 3/5 - sigma), sigma = {float(s):.0e}: r = {mp.nstr(ratio([s, Fr(4, 10), 1 - s - Fr(4, 10)]), 10)}')
random.seed(3)
worst = (0, None)
for tr in range(300):
    k = random.randint(3, 7)
    m = tr % 5
    if m == 0:
        raw = [Fr(random.randint(1, 50)) for _ in range(k)]
    elif m == 1:
        raw = [Fr(1, 3 ** random.randint(0, 8)) for _ in range(k)]
    elif m == 2:
        raw = [Fr(random.randint(1, 50)) for _ in range(k)]; raw[0] = Fr(1, 10 ** random.randint(2, 6))
    elif m == 3:
        raw = [Fr(random.randint(1, 50)) for _ in range(k)]; raw[0] = Fr(10 ** random.randint(2, 6))
    else:
        raw = [Fr(1, 10 ** random.randint(2, 6)) for _ in range(k)]
        raw[random.randrange(1, k)] = Fr(1); raw[random.randrange(1, k)] = Fr(1)
    tot = sum(raw)
    xs = [v / tot for v in raw]
    r = ratio(xs)
    if r > worst[0]:
        worst = (r, [float(v) for v in xs])
print(f'largest r over 300 rational weight vectors: {mp.nstr(worst[0], 10)} at x = {[round(v, 6) for v in worst[1]]}')
print('numerical exploration done')
