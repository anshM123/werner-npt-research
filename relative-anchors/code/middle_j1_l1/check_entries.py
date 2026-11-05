"""The closed forms of l1sym.Model, entry by entry, against the letter-1 block of (R_inf)_1 built from Lemma B by
j1lib.py (rescaled frame, units 2^{k-1}), exactly, at rational two-level weights, with the B-part scaled by 1 and by
1001/1000:
  PP entries Y~[(r1,c1),(r2,c2)] (every pair of P-states);
  PQ entries Y~[(r,c),{a,b}] (every P-state and pair);
  the vacuum block: Y~_QQ = 2 W^{-1} for b = 1, and 2 D2 + b (2 W^{-1} - 2 D2) in general, with W built from Model.W.
Also the Lemma G facts used: W is symmetric and positive definite (exact LDL), so the vacuum block is positive
definite.  Here the copies are concrete (labels (level, None, copy)), so every equality pattern is exercised.
Usage: python check_entries.py"""
import time, itertools
from fractions import Fraction as Fr
from j1lib import J1, block_states, mu
from l1sym import Model, minv

t0 = time.time()


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


ok = True
pts = [(Fr(1), Fr(1, 3), 2, 2), (Fr(1), Fr(5, 2), 1, 3), (Fr(3), Fr(1), 3, 1), (Fr(1, 2), Fr(3), 3, 2),
       (Fr(1), Fr(1, 7), 2, 4), (Fr(2), Fr(9, 4), 4, 2)]
for (w0, w1, u0, u1) in pts:
    for b in (Fr(1), Fr(1001, 1000)):
        k = u0 + u1
        w = [w0] * u0 + [w1] * u1
        lab = [(0 if t < u0 else 1, None, t) for t in range(k)]
        J = J1(w, bscale=b)
        st = block_states(k, ())
        Mb = J.matrix(st)
        m = [mu(z, J.w) for z in st]
        Mt = [[m[i] * Mb[i][j] * m[j] for j in range(len(st))] for i in range(len(st))]
        Mo = Model(Fr(1), [w0 * w0, w1 * w1], [Fr(u0), Fr(u1)], b)
        Pi = [i for i, (a, y) in enumerate(st) if a is not None]
        Qi = [i for i, (a, y) in enumerate(st) if a is None]
        rc = {i: (st[i][0][1], st[i][1].index(1)) for i in Pi}
        pr = {i: tuple(t for t in range(k) if st[i][1][t] == 1) for i in Qi}
        good_pp = all(Mt[i][j] == Mo.PP(lab[rc[i][0]], lab[rc[i][1]], lab[rc[j][0]], lab[rc[j][1]])
                      for i in Pi for j in Pi)
        good_pq = all(Mt[i][j] == Mo.PQ(lab[rc[i][0]], lab[rc[i][1]], lab[pr[j][0]], lab[pr[j][1]])
                      for i in Pi for j in Qi)
        Wm = [[Mo.W(lab[pr[i][0]], lab[pr[i][1]], lab[pr[j][0]], lab[pr[j][1]]) for j in Qi] for i in Qi]
        D2 = [[Mo.D2(lab[pr[i][0]], lab[pr[i][1]], lab[pr[j][0]], lab[pr[j][1]]) for j in Qi] for i in Qi]
        Wi = minv(Wm, Fr(1))
        QQ = [[2 * D2[a][c] + b * (2 * Wi[a][c] - 2 * D2[a][c]) for c in range(len(Qi))] for a in range(len(Qi))]
        good_qq = all(Mt[Qi[a]][Qi[c]] == QQ[a][c] for a in range(len(Qi)) for c in range(len(Qi)))
        wsym = all(Wm[a][c] == Wm[c][a] for a in range(len(Qi)) for c in range(len(Qi)))
        wpd = ldl_pd(Wm)
        good = good_pp and good_pq and good_qq and wsym and wpd
        ok &= good
        print(f'w = ({w0}, {w1}), (u0, u1) = ({u0}, {u1}), b = {b}: PP {good_pp}, PQ {good_pq}, QQ {good_qq}, '
              f'W symmetric {wsym}, W positive definite {wpd} [{time.time() - t0:.1f}s]', flush=True)
print('PASS' if ok else 'FAIL')
