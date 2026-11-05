"""Direct exact checks of the conclusion on the whole letter-1 block of (R_inf)_1 (j1lib.py, Lemma B), at two-level
rational weights, and agreement with the reduced forms.

For k = 4..7, splits (u0, u1) with u0, u1 >= 1 and weights w0 = 1 on level 0, w1 on level 1:
 (i)  the indicators 1_c of the column groups (the kernel vectors p (x) |1@c> in the rescaled frame) are annihilated
      by the rescaled block;
 (ii) the gauge-fixed block (the k diagonal states P_{c,c} deleted; each 1_c has coefficient 1 on P_{c,c} and 0 on the
      other deleted states) is positive definite (exact LDL), so the block is >= 0 with kernel exactly span{1_c};
 (iii) the reduced forms of l1sym.py agree: every part that exists at (u0, u1) (irreps on either level) has a
      positive definite vacuum block and a positive definite gauge-fixed Schur complement.
These are exact checks at points, not proofs; the proof is the certificate (cert.py).
Usage: python check_block.py"""
import time
from fractions import Fraction as Fr
from j1lib import J1, block_states, mu
from l1sym import Model, reduce_part, keep_for, det

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


def block_verdict(w, b=Fr(1)):
    k = len(w)
    J = J1(w, bscale=b)
    st = block_states(k, ())
    M = J.matrix(st)
    m = [mu(z, J.w) for z in st]
    Mt = [[m[i] * M[i][j] * m[j] for j in range(len(st))] for i in range(len(st))]
    kerok = True
    for c in range(k):
        v = [Fr(1) if (a is not None and y[c] == 1) else Fr(0) for (a, y) in st]
        kerok &= all(sum(Mt[i][j] * v[j] for j in range(len(st)) if v[j]) == 0 for i in range(len(st)))
    keep = [i for i, (a, y) in enumerate(st) if not (a is not None and y[a[1]] == 1)]
    pd = ldl_pd([[Mt[i][j] for j in keep] for i in keep])
    return kerok, pd


def parts_at(u0, u1):
    """(part, levels exchanged?) for every part that exists at (u0, u1)."""
    out = [('inv', False)]
    for sw, (a, c) in ((False, (u0, u1)), (True, (u1, u0))):
        if a >= 2:
            out.append(('std', sw))
        if a >= 4:
            out.append(('two', sw))
        if a >= 3:
            out.append(('anti', sw))
    if u0 >= 2 and u1 >= 2:
        out.append(('mix', False))
    return out


def reduced_verdict(w0, w1, u0, u1, b=Fr(1)):
    ok = True
    for part, sw in parts_at(u0, u1):
        a, c = (u1, u0) if sw else (u0, u1)
        qa, qc = (w1 * w1, w0 * w0) if sw else (w0 * w0, w1 * w1)
        kp, kq, gauge = keep_for(part, a, c)
        red = reduce_part(Model(Fr(1), [qa, qc], [Fr(a), Fr(c)], b), part, kp, kq)
        S = red['S']
        keep = [i for i in range(len(S)) if i not in gauge]
        Sg = [[S[i][j] for j in keep] for i in keep]
        ok &= ldl_pd(Sg) and (ldl_pd(red['RQQ']) if 'RQQ' in red else True)
    return ok


if __name__ == '__main__':
    ok = True
    cnt = 0
    cases = []
    for k in (4, 5, 6):
        for u0 in range(1, k):
            for w1 in (Fr(1, 40), Fr(1, 3), Fr(7, 2), Fr(30)):
                cases.append((u0, k - u0, w1))
    for (u0, u1) in ((4, 3), (2, 5), (6, 1)):
        for w1 in (Fr(1, 20), Fr(5)):
            cases.append((u0, u1, w1))
    for (u0, u1, w1) in cases:
        w = [Fr(1)] * u0 + [w1] * u1
        kerok, pd = block_verdict(w)
        rv = reduced_verdict(Fr(1), w1, u0, u1)
        good = kerok and pd and rv
        ok &= good
        cnt += 1
        print(f'k = {u0 + u1}, (u0, u1) = ({u0}, {u1}), w1 = {w1}: kernel vectors {kerok}, gauge-fixed block PD {pd}, '
              f'reduced forms PD {rv} [{time.time() - t0:.1f}s]', flush=True)
    print(f'{cnt} weight vectors; time {time.time() - t0:.1f}s')
    print('PASS' if ok else 'FAIL')
