"""Proof of (B_2) for u = 2: lambda_min(Ghat_tt) > H (1 - s_t/sigma), H = 2/5, for every m >= 1 and all weights.

Setting (build_B.py): N = G_z - h diag(1, s_t (diag(s_U) - s_U s_U^T/(sigma + S))), h = H (1 - th), variables
(th, mu, la) in D = (0,1] x (0,1) x (0,1), N = Nt/den.  Ghat_tt - h I > 0 iff N > 0.  Steps:
  (1) den > 0 on D: den is a product of factors; each factor is certified > 0 on the closed cube by Bernstein
      coefficients, or is a power of a variable factor that is positive on D.
  (2) det N = m3/den^3 with m3 = th^2 mu^2 la (1-la) (1-mu)^6 P3; P3 > 0 on [0,1]^3 minus the edge E1 = {th = 1,
      mu = 0}: two blow-up charts (th, mu) = (1 - e, e v) and (1 - mu v, mu); in each the substituted P3 is divided by
      the largest power of e (resp. mu) and the quotient is certified > 0 on [0,1]^3 (tensor Bernstein, exact).
  (3) N > 0 at one point of D (exact LDL^T).
  (4) D is connected and det N != 0 on D, so the inertia of N is constant on D; by (3), N > 0 on D."""
import pickle, sys, time, flint
from math import lcm
from fractions import Fraction as Fr
from bern import from_dict, certify
import build_B as bb

H = Fr(sys.argv[1]) if len(sys.argv) > 1 else Fr(2, 5)
th, mu, la = bb.th, bb.mu, bb.la
ok = True
t00 = time.time()


def cert(name, p, gens_ctx_degs=None, corner=None):
    dct = p.to_dict()
    L = 1
    for c in dct.values():
        L = lcm(L, int(c.q))
    idct = {tuple(int(e) for e in m): int(c * L) for m, c in dct.items()}
    degs = tuple(int(d) for d in p.degrees())
    t0 = time.time()
    good, st = certify(from_dict(idct, degs), max_boxes=20000, corner=corner)
    print(f'  {name}: degrees {degs}, {len(dct)} terms: {"certified > 0 on [0,1]^3" if good else "FAILED"} {st} ({time.time() - t0:.1f} s)', flush=True)
    return good


# ---- (1) the denominator
num, den = bb.load()
c, fac = den.factor()
print(f'(1) den = {c} * product of {len(fac)} factors', flush=True)
sign = 1 if c > 0 else -1
for f_, e_ in fac:
    # each factor: certify f > 0 or -f > 0 on the closed cube (then it has constant sign on D)
    if cert(f'factor (degrees {f_.degrees()}) ^ {e_}', f_):
        s_ = 1
    elif cert(f'minus factor (degrees {f_.degrees()}) ^ {e_}', -f_):
        s_ = -1
    elif f_ in (th, mu, la, 1 - la, 1 - mu, 1 - th, mu - 1, la - 1, th - 1):
        s_ = -1 if f_ in (mu - 1, la - 1, th - 1) else 1
        print(f'  factor {f_}: sign {s_} on D (linear, vanishes only on the boundary)')
    else:
        s_ = 0
    ok &= s_ != 0
    sign *= s_ ** e_
print(f'  den has sign {sign} on D')
ok &= sign == 1

# ---- (2) det N
Hq = flint.fmpq(H.numerator, H.denominator)
out, _ = bb.build(H)
name, m3, rem = out[2]
print(f'(2) det: m3 = {rem} * P3, P3 degrees {m3.degrees()}', flush=True)
ok &= all(nm in ('th', 'mu', 'la', '1-la', '1-mu') for nm, e in rem)
ctx3 = bb.ctx
e_, v_, l_ = ctx3.gens()      # reuse the same context: (th, mu, la) slots now hold (e, v, la) and (mu, v, la)
for chart, sub, div in (('alpha: th = 1 - e, mu = e v', (1 - e_, e_ * v_, l_), e_),
                        ('beta: th = 1 - mu v, mu = mu', (1 - e_ * v_, e_, l_), e_)):
    q = m3.compose(*sub, ctx=ctx3)
    k = 0
    while True:
        qq, r = divmod(q, div)
        if r != 0:
            break
        q = qq; k += 1
    print(f'  chart {chart}: divided by (first variable)^{k}', flush=True)
    ok &= cert(f'P3 chart {chart.split(":")[0]}', q)

# ---- (3) one point
pt = (Fr(1, 2), Fr(1, 2), Fr(1, 3))
q = [flint.fmpq(x.numerator, x.denominator) for x in pt]
sig = 1 - pt[1]; st = pt[0] * sig; sU = [pt[1] * pt[2], pt[1] * (1 - pt[2])]; h = H * (1 - pt[0])
Gv = [[Fr(int((num[i][j](*q) / den(*q)).p), int((num[i][j](*q) / den(*q)).q)) for j in range(3)] for i in range(3)]
Oinv = [[Fr(1), Fr(0), Fr(0)], [Fr(0), st * (sU[0] - sU[0] ** 2), -st * sU[0] * sU[1]], [Fr(0), -st * sU[0] * sU[1], st * (sU[1] - sU[1] ** 2)]]
Nv = [[Gv[i][j] - h * Oinv[i][j] for j in range(3)] for i in range(3)]
piv = []
A = [row[:] for row in Nv]
for p in range(3):
    piv.append(A[p][p])
    for i in range(p + 1, 3):
        f = A[i][p] / A[p][p]
        for j in range(p + 1, 3):
            A[i][j] -= f * A[p][j]
print(f'(3) N at (th, mu, la) = {tuple(str(x) for x in pt)}: LDL pivots {[float(x) for x in piv]}', flush=True)
ok &= all(x > 0 for x in piv)
print(f'total {time.time() - t00:.1f} s')
print('PASS' if ok else 'FAIL')
