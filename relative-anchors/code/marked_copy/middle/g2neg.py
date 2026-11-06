"""Negative control for the determinant route: with B replaced by (1 + 1/1000) B, (a) the kernel vectors are no
longer annihilated (so the gauge argument is not available), and (b) at the slice (t, i) = (4, 7) the determinant of
the gauge-fixed form over Q(q0, s), after q0 = 4 + x, has coefficients of both signs and takes a negative value at
an explicit point, so the certificate of g2final.py cannot pass.  Expected output: FAIL (the control succeeds)."""
import time
import flint
import famlib
famlib.SMIN['q0'] = 4
famlib.SMIN['t'] = 4
famlib.SMIN['i'] = 7
import splitlib
from fractions import Fraction as Fr
from prove_fam import fam_kernel
from flintfield import Field, F

t0 = time.time()
C = famlib.FamCalc(None, None, None, bscale=Fr(1001, 1000))
names, Q = famlib.reduced(C)
vecs, split = splitlib.basis(names, C.q0, C.m, C.i, C.one)
two, three, M = splitlib.congruence(names, Q, vecs, C.zero)
fk = fam_kernel(names, C)
annihilated = True
for lab, v, dl in fk:
    w = splitlib.to_new(v, split)
    for a in two:
        s = C.zero
        for nm, cf in w.items():
            s = s + M[a, nm] * cf
        if s != 0:
            annihilated = False
print(f'(a) kernel vectors annihilated with B scaled: {annihilated}', flush=True)
keep = [nm for nm in two if nm not in ('U10', 'U11')]
K = Field(['q0', 's'])
q, s = K.ctx.gens()
c4, c7 = K.ctx.from_dict({(0, 0): 4}), K.ctx.from_dict({(0, 0): 7})
sp = lambda P: P.compose(q, c4, c7, s, ctx=K.ctx)
A = [[(F(K, sp(M[a, b].n), sp(M[a, b].d)) if M[a, b] != 0 else K(0)) for b in keep] for a in keep]
n = len(A)
det = K(1)
for c in range(n):
    p = A[c][c]
    det = det * p
    for r in range(c + 1, n):
        if A[r][c] != 0:
            f = A[r][c] / p
            for cc in range(c + 1, n):
                A[r][cc] = A[r][cc] - f * A[c][cc]
x = K.ctx.gens()[0]
Nsh = det.n.compose(x + 4, s, ctx=K.ctx)
co = list(Nsh.to_dict().values())
npos, nneg = sum(1 for c in co if c > 0), sum(1 for c in co if c < 0)
print(f'(b) det at (t, i) = (4, 7) with B scaled, after q0 = 4 + x: {npos} positive, {nneg} negative coefficients',
      flush=True)
neg_pt = None
for qv in (4, 5, 8, 20, 100):
    for sv in (Fr(1, 10 ** 6), Fr(1, 1000), Fr(1, 10), Fr(1), Fr(10), Fr(1000), Fr(10 ** 6)):
        val = det.n(flint.fmpq(qv), flint.fmpq(sv.numerator, sv.denominator)) / \
            det.d(flint.fmpq(qv), flint.fmpq(sv.numerator, sv.denominator))
        if val < 0 and neg_pt is None:
            neg_pt = (qv, sv)
print(f'    a point (q0, s) with negative determinant: {neg_pt} [{time.time()-t0:.1f}s]')
passes = annihilated and nneg == 0 and neg_pt is None
print('PASS' if passes else 'FAIL')
