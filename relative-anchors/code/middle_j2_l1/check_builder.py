"""The builder l2block.py against the Lemma B builder j1lib.py (j = 2, units 2^{k-2}), exact, entry by entry.

At rational weight vectors with k = 5, 6, 7 (general weights, and two values with the levels interleaved), with the
B-part scaled by 1 and by 1001/1000: the rescaled j1lib block (entries times mu(row) mu(column)) equals l2block.full().
Also, at b = 1: the vacuum identity Y_QQ (What3 / 2) = I (Lemma G with Woodbury), homogeneity of degree 3 in the
squared weights (Y(c s) = c^3 Y(s)), the kernel indicators 1_X are in the kernel, the nullity is C(k, 2), and the block
with the states (min X, X) deleted is positive definite (exact LDL).  Usage: python check_builder.py"""
import time
from fractions import Fraction as Fr
from math import comb
from j1lib import J1, mu
import flint
from r2lib import states_L1, ldl_pd
from l2block import Block

t0 = time.time()
ok = True
cases = [
    [Fr(1), Fr(2), Fr(3), Fr(1, 2), Fr(5, 3)],
    [Fr(1), Fr(1, 3), Fr(1), Fr(1, 3), Fr(1)],
    [Fr(2), Fr(2), Fr(1), Fr(1), Fr(2)],
    [Fr(1, 2), Fr(3), Fr(1, 2), Fr(3), Fr(3), Fr(1, 2)],
    [Fr(1), Fr(1), Fr(1), Fr(1), Fr(1), Fr(7)],
    [Fr(3), Fr(1, 5), Fr(2), Fr(1), Fr(4, 3), Fr(5, 2)],
    [Fr(1), Fr(1, 3), Fr(1), Fr(1, 3), Fr(1), Fr(1, 3), Fr(1)],
    [Fr(2), Fr(5), Fr(5), Fr(2), Fr(5), Fr(2), Fr(2)],
]
for w in cases:
    k = len(w)
    s = [x * x for x in w]
    st = states_L1(k)
    for b in (Fr(1), Fr(1001, 1000)):
        J = J1(w, bscale=b)
        M = J.matrix(st)
        mus = [mu(x, w) for x in st]
        B = Block(s, Fr(1), bscale=b)
        F = B.full()
        n = len(st)
        eq = all(F[i][j] == M[i][j] * mus[i] * mus[j] for i in range(n) for j in range(n))
        ok &= eq
        line = f'k={k} w={[str(x) for x in w]} b={b}: l2block = rescaled j1lib: {eq}'
        if b == 1:
            QQ, W3 = B.QQ(), B.What3()
            m = len(QQ)
            wood = all(sum(QQ[i][l] * W3[l][j] for l in range(m)) / 2 == (1 if i == j else 0)
                       for i in range(m) for j in range(m))
            c = Fr(7, 3)
            F2 = Block([c * x for x in s], Fr(1)).full()
            hom = all(F2[i][j] == c ** 3 * F[i][j] for i in range(n) for j in range(n))
            nP = len(B.P)
            inker = True
            for X in B.sets[2]:
                v = [1 if (i < nP and B.P[i][1] == X) else 0 for i in range(n)]
                inker &= all(sum(F[i][j] * v[j] for j in range(n)) == 0 for i in range(n))
            null = n - flint.fmpq_mat(n, n, [flint.fmpq(x.numerator, x.denominator) for row in F for x in row]).rank()
            keep = [i for i in range(n) if not (i < nP and B.P[i][0] == B.P[i][1][0])]
            pd = ldl_pd([[F[i][j] for j in keep] for i in keep])
            good = wood and hom and inker and null == comb(k, 2) and pd
            ok &= good
            line += (f'; Y_QQ What3/2 = I: {wood}; homogeneous of degree 3: {hom}; 1_X in kernel: {inker}; '
                     f'nullity {null} (C(k,2) = {comb(k, 2)}); gauge-fixed PD: {pd}')
        print(line + f' [{time.time() - t0:.1f}s]', flush=True)
print('PASS' if ok else 'FAIL')
