"""The Lemma B builder j1lib.py agrees exactly with ylib.Builder (the exact builder of relative-anchors/code, rescaled
basis) on the type I blocks T = {t} of (R_inf)_1, at rational weights (equal U-weights, two U-weight levels, generic),
k = 3, ..., 6.  Also: p (x) |2@t> (the vector w on the P-states) is in the kernel.  Usage: python check_j1_ylib.py"""
import time
from fractions import Fraction as Fr
import flint
from ylib import Builder, point_ring
from j1lib import J1, block_states, mu, to_ylib_state

t0 = time.time()
ok = True
cases = [([Fr(2), Fr(1), Fr(1)], 0), ([Fr(1, 3), Fr(2), Fr(2), Fr(2)], 0), ([Fr(3), Fr(1, 2), Fr(1, 2), Fr(5, 2), Fr(5, 2)], 0),
         ([Fr(1), Fr(2), Fr(3), Fr(5), Fr(1, 3)], 2), ([Fr(7, 5), Fr(1), Fr(1), Fr(1), Fr(4), Fr(4)], 0),
         ([Fr(3, 2), Fr(1), Fr(1, 4), Fr(2), Fr(5, 2), Fr(2, 3)], 4)]
for w, t in cases:
    k = len(w)
    Bd = Builder(point_ring([flint.fmpq(x.numerator ** 2, x.denominator ** 2) for x in w]))
    J = J1(w)
    st = block_states(k, (t,))
    M = J.matrix(st)
    ys = [to_ylib_state(x, k) for x in st]
    good = True
    for i in range(len(st)):
        for j in range(len(st)):
            yv = Bd.Y(ys[i], ys[j])
            mine = M[i][j] * mu(st[i], w) * mu(st[j], w) / 2 ** (k - 1)
            if Fr(int(yv.p), int(yv.q)) != mine:
                good = False
    v = [w[a[1]] if (a is not None and a[0] == 1) else 0 for (a, y) in st]
    kern = all(sum(M[i][j] * v[j] for j in range(len(st))) == 0 for i in range(len(st)))
    print(f'k = {k}, t = {t}, w = {[str(x) for x in w]}: size {len(st)} (3k - 1 = {3 * k - 1}), j1lib = ylib exactly: '
          f'{good}, p (x) |2@t> in the kernel: {kern} [{time.time() - t0:.1f}s]', flush=True)
    ok &= good and kern and len(st) == 3 * k - 1
print('PASS' if ok else 'FAIL')
