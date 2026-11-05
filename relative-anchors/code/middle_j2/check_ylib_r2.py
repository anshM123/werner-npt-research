"""The Lemma B builder j1lib.py (used with j = 2, units 2^{k-2}) equals the builder from the definitions ylib.py
(M, V_p, D_p, R_p and the partial transpose, rescaled basis) exactly, entry by entry, on the blocks of the three
families of (R_inf)_2: k = 5 at two rational weight vectors (all three families, also with two weight values), and
k = 6 at one weight vector (families 1 and 2).  Relation: Y'(ylib) = mu_X mu_X' 2^{-(k-2)} M(j1lib).
Usage: python check_ylib_r2.py"""
import time
from fractions import Fraction as Fr
import flint
from ylib import Builder, point_ring
from j1lib import J1, mu, to_ylib_state
from r2lib import states_T2, states_T1, states_L1

t0 = time.time()
ok = True
cases = [(5, [Fr(1), Fr(2), Fr(3), Fr(5), Fr(1, 3)], (1, 2, 3)),
         (5, [Fr(1, 2), Fr(3), Fr(3), Fr(1, 2), Fr(1, 2)], (1, 2, 3)),
         (6, [Fr(2), Fr(1), Fr(1, 3), Fr(1), Fr(1, 3), Fr(1)], (1, 2))]
for k, w, fams in cases:
    Bd = Builder(point_ring([flint.fmpq(x.numerator ** 2, x.denominator ** 2) for x in w]))
    J = J1(w)
    Z = k - 2
    for fam in fams:
        st = states_T2(k, 0, 1) if fam == 1 else (states_T1(k, 0) if fam == 2 else states_L1(k))
        M = J.matrix(st)
        ys = [to_ylib_state(x, k) for x in st]
        good = True
        for i in range(len(st)):
            for j in range(i, len(st)):
                yv = Bd.Y(ys[i], ys[j])
                mine = M[i][j] * mu(st[i], w) * mu(st[j], w) / 2 ** Z
                if Fr(int(yv.p), int(yv.q)) != mine:
                    good = False
        ok &= good
        print(f'k={k} w={[str(x) for x in w]} family {fam}: {len(st)} states, j1lib = ylib exactly: {good} '
              f'[{time.time() - t0:.1f}s]', flush=True)
print('PASS' if ok else 'FAIL')
