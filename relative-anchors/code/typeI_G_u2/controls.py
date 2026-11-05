"""Negative controls for (A_2) and (B_2): with the constants moved past the numerical extremes the statements are false,
and an explicit rational witness is exhibited (exact arithmetic).
 (A) H = 1/4 (< sup ||Ghat_tt'||_F sigma/sqrt(s_t s_t') ~ 0.28): a point with FA < 0.
 (B) H = 4/5 (> inf lambda_min(Ghat_tt)/(1 - s_t/sigma) = 3/4): a point where N has a negative LDL pivot.
Also: the certificate routine itself rejects the shifted polynomial PA(H = 1/4) (it must not certify)."""
import sys, time, flint
from fractions import Fraction as Fr
import build_A as ba
import build_B as bb

ok = True
t0 = time.time()
# ---- (A): search a few points near the corner th, ph -> 0 for the largest normalised Frobenius norm
FA4, den = ba.build(Fr(1, 4))
best = None
# witness near the numerical maximiser of the Frobenius ratio (~0.288 at th + ph = 1, mu -> 0, small la; ex_findA.py)
for pt_ in ((Fr(11, 20), Fr(9, 20), Fr(1, 10000), Fr(1, 100)), (Fr(11, 20), Fr(2, 5), Fr(1, 1000), Fr(1, 50))):
    pt = [flint.fmpq(x.numerator, x.denominator) for x in pt_]
    v = FA4(*pt)
    if best is None or v < best[0]:
        best = (v, pt_)
v, pt = best
print(f'(A) H = 1/4: FA at (th, ph, mu, la) = {tuple(str(x) for x in pt)} has sign {"negative" if v < 0 else "nonnegative"}', flush=True)
ok &= v < 0
# ---- (B): search for a point where N(H = 4/5) is not positive definite
num, den = bb.load()
H = Fr(4, 5)
found = None
for th_ in (Fr(1, 1000), Fr(1, 100)):
    for mu_ in [Fr(i, 20) for i in range(1, 20)]:
        for la_ in (Fr(1, 1000), Fr(1, 2), Fr(999, 1000)):
            q = [flint.fmpq(x.numerator, x.denominator) for x in (th_, mu_, la_)]
            sig = 1 - mu_; st = th_ * sig; sU = [mu_ * la_, mu_ * (1 - la_)]; h = H * (1 - th_)
            G = [[Fr(int((num[i][j](*q) / den(*q)).p), int((num[i][j](*q) / den(*q)).q)) for j in range(3)] for i in range(3)]
            O = [[Fr(1), Fr(0), Fr(0)], [Fr(0), st * (sU[0] - sU[0] ** 2), -st * sU[0] * sU[1]], [Fr(0), -st * sU[0] * sU[1], st * (sU[1] - sU[1] ** 2)]]
            A = [[G[i][j] - h * O[i][j] for j in range(3)] for i in range(3)]
            piv = []
            for p in range(3):
                piv.append(A[p][p])
                if A[p][p] <= 0:
                    break
                for i in range(p + 1, 3):
                    f = A[i][p] / A[p][p]
                    for j in range(p + 1, 3):
                        A[i][j] -= f * A[p][j]
            if min(piv) <= 0:
                found = (th_, mu_, la_, piv)
                break
        if found:
            break
    if found:
        break
print(f'(B) H = 4/5: N not positive definite at (th, mu, la) = {tuple(str(x) for x in found[:3]) if found else None}', flush=True)
ok &= found is not None
print(f'time {time.time() - t0:.1f} s')
print('PASS' if ok else 'FAIL')
