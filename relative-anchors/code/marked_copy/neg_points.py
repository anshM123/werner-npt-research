"""Negative controls at explicit points, from the Lemma B builder alone.
For k = 3..KMAX, every j = -1..k and i = 0: at sigma = 10^-6 (10^6 for j = -1), both rational squares, take the
gauge of the certificate (cert_fixed.choose_gauge), the kept explicit twisted vectors v_1..v_n (explicit.py), and the
rescaled block built by j1lib.py with the B-part scaled by 1 and by 1001/1000.  The restriction V^T Y~ V must be
positive definite for the true block (exact LDL) and not positive definite for the scaled block.
Usage: python neg_points.py KMAX"""
import sys, time
from fractions import Fraction as Fr
import flint
import mlib, explicit
from cert_fixed import choose_gauge
from ufield import RF

KMAX = int(sys.argv[1])
T0 = time.time()
ok = True
fq = flint.fmpq


def pd(M):
    """exact: True iff the symmetric rational matrix M is positive definite (all LDL pivots > 0)"""
    n = len(M)
    A = [[fq(x.numerator, x.denominator) for x in row] for row in M]
    for t in range(n):
        piv = A[t][t]
        if piv <= 0:
            return False
        for a in range(t + 1, n):
            f = A[a][t] / piv
            if f != 0:
                for b in range(t + 1, n):
                    A[a][b] -= f * A[t][b]
    return True


cnt = 0
for k in range(3, KMAX + 1):
    for j in range(-1, k + 1):
        w0 = Fr(1000) if j == -1 else Fr(1, 1000)
        sig = w0 * w0
        w = [w0] + [Fr(1)] * (k - 1)
        i = 0
        names, norms, Q = mlib.reduced(mlib.Calc(k, j, i, sigma=None))
        if not names:
            continue
        kv = mlib.kernel_vectors(names, k, j, i, RF(mlib.PL([0, 1])))
        kvR = [(lab, {nm: (cf if isinstance(cf, RF) else RF(fq(cf))) for nm, cf in v.items()}) for lab, v in kv]
        S, how = choose_gauge(names, kvR)
        keep = [names[c] for c in range(len(names)) if c not in S]
        st, vecs = explicit.twisted(k, j, i)
        res = []
        for bs in (Fr(1), Fr(1001, 1000)):
            st, M = explicit.block(k, j, w, bscale=bs)
            R = [[explicit.form(M, vecs[a], vecs[b]) for b in keep] for a in keep]
            res.append(pd(R))
        cnt += 1
        if res != [True, False]:
            ok = False
            print(f'k={k} j={j} sigma={sig}: true block PD on the complement: {res[0]}, scaled block PD: {res[1]}')
    print(f'k={k}: done ({time.time() - T0:.0f}s)', flush=True)
print(f'{cnt} points: the true restricted block is positive definite and the scaled one is not')
print('PASS' if ok else 'FAIL', f'{time.time() - T0:.0f}s')
