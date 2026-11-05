"""Step 1: the direct Lemma-B implementation j1lib.py agrees exactly with ylib.Builder (rescaled basis) on every block
of (R_inf)_1 that is not of type II (letter-1, T = {t}, T = {t, t'}), at rational weights (w rational, so s = w^2),
k = 4, 5 (all three blocks) and k = 6 (all three blocks, one weight vector).  Also: block sizes, the kernel vectors
p (x) y, and that the block is PSD with kernel exactly p (x) (columns of the block) (exact symmetric elimination in a
gauge).  Usage: python check_ylib.py"""
import sys, time
from fractions import Fraction as Fr
import flint
from ylib import Builder, point_ring
from j1lib import J1, block_states, mu, to_ylib_state

t0 = time.time()
ok = True


def gauge_pd(M, states, F=Fr):
    """delete (1@g, y) for g = 1 (one per column y that carries rows (1@r, y)); exact LDL of the rest; return
    (all pivots > 0, number of deleted)."""
    keep = [i for i, (a, y) in enumerate(states) if not (a is not None and a[0] == 1 and a[1] == 1)]
    A = [[M[i][j] for j in keep] for i in keep]
    n = len(A)
    for c in range(n):
        p = A[c][c]
        if p <= 0:
            return False, len(states) - len(keep)
        for i in range(c + 1, n):
            if A[i][c] != 0:
                f = A[i][c] / p
                for j in range(c + 1, n):
                    A[i][j] -= f * A[c][j]
    return True, len(states) - len(keep)


cases = [(4, [Fr(1), Fr(2), Fr(3), Fr(5)]), (4, [Fr(7, 3), Fr(1, 2), Fr(1), Fr(4)]),
         (5, [Fr(1), Fr(2), Fr(3), Fr(5), Fr(1, 3)]), (5, [Fr(3, 2), Fr(1), Fr(1, 4), Fr(2), Fr(5, 2)]),
         (6, [Fr(1), Fr(2), Fr(1, 2), Fr(3), Fr(5, 3), Fr(2, 3)])]
for k, w in cases:
    Bd = Builder(point_ring([flint.fmpq(x.numerator ** 2, x.denominator ** 2) for x in w]))
    J = J1(w)
    Z = k - 1
    for T in [(), (0,), (0, 1), (2,), (1, 3)]:
        st = block_states(k, T)
        M = J.matrix(st)
        # compare with ylib: Y' = mu Y mu, Y = 2^{-Z} M
        ys = [to_ylib_state(x, k) for x in st]
        maxdiff = 0
        good = True
        for i in range(len(st)):
            for j in range(len(st)):
                yv = Bd.Y(ys[i], ys[j])
                mine = M[i][j] * mu(st[i], w) * mu(st[j], w) / 2 ** Z
                if Fr(int(yv.p), int(yv.q)) != mine:
                    good = False
        # kernel p (x) y
        cols = sorted({y for (a, y) in st if a is not None and a[0] == 1})
        kern_ok = True
        for y in cols:
            v = [w[a[1]] if (a is not None and a[0] == 1 and yy == y) else 0 for (a, yy) in st]
            Mv = [sum(M[i][j] * v[j] for j in range(len(st))) for i in range(len(st))]
            kern_ok &= all(x == 0 for x in Mv)
        pd, ndel = gauge_pd(M, st)
        expect = {0: (3 * k * k - k) // 2, 1: 3 * k - 1, 2: 3}[len(T)]
        print(f"k={k} w={[str(x) for x in w]} T={T}: size {len(st)} (expected {expect}), equal to ylib exactly: {good}, "
              f"p(x)y in kernel ({len(cols)} vectors): {kern_ok}, gauge-fixed block PD: {pd} [{time.time()-t0:.1f}s]",
              flush=True)
        ok &= good and kern_ok and pd and len(st) == expect and ndel == len(cols)
print("PASS" if ok else "FAIL")
