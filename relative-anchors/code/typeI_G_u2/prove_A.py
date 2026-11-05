"""Proof of (A_2) for u = 2: ||Ghat_tt'||_2 <= ||Ghat_tt'||_F <= H sqrt(s_t s_t')/sigma, H = 2/5, all weights.

FA (build_A.py) = H^2 s_t^2 s_t'^2 (s_3 s_4)^2 den^2 - tr(Om~_t Nm Om~_t' Nm^T) in (th, ph, mu, la), th + ph <= 1.
FA = th^a ph^b ... * PA (powers of factors positive on the open domain removed).  The triangle th + ph <= 1 is
mapped onto the square by th = x, ph = (1 - x) y (onto); PA(x, (1-x) y, mu, la) is certified > 0 on [0,1]^4 by
exact tensor Bernstein coefficients.  Hence FA > 0, i.e. ||Ghat_tt'||_F^2 < H^2 s_t s_t'/sigma^2, on the whole
domain (den != 0 there: den is the same product of factors as in prove_B plus factors certified below)."""
import pickle, sys, time, flint
from math import lcm
from fractions import Fraction as Fr
from bern import from_dict, certify
import build_A as ba

H = Fr(sys.argv[1]) if len(sys.argv) > 1 else Fr(2, 5)
max_boxes = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
th, ph, mu, la = ba.th, ba.ph, ba.mu, ba.la
ok = True
t00 = time.time()


def cert(name, p, mb=max_boxes):
    dct = p.to_dict()
    L = 1
    for c in dct.values():
        L = lcm(L, int(c.q))
    idct = {tuple(int(e) for e in m): int(c * L) for m, c in dct.items()}
    degs = tuple(int(d) for d in p.degrees())
    t0 = time.time()
    good, st = certify(from_dict(idct, degs), max_boxes=mb)
    print(f'  {name}: degrees {degs}, {len(dct)} terms: {"certified > 0 on [0,1]^4" if good else "FAILED"} {st} ({time.time() - t0:.1f} s)', flush=True)
    return good


FA, den = ba.build(H)
# the denominator of G01: every factor nonzero on the open domain (certified on the closed cube after the triangle map)
c, fac = den.factor()
print(f'den = {c} * {len(fac)} factors', flush=True)
x, y, m_, l_ = ba.ctx.gens()
for f_, e_ in fac:
    g = f_.compose(x, (1 - x) * y, m_, l_, ctx=ba.ctx)
    if f_ == mu - 1 or f_ == 1 - mu:
        print(f'  factor {f_} ^ {e_}: nonzero for mu < 1 (sigma > 0)')
        continue
    good = cert(f'den factor ^{e_} (deg {f_.total_degree()})', g, 200) or cert(f'minus den factor ^{e_}', -g, 200)
    ok &= good
q, rem = ba.strip(FA)
print(f'FA = {rem} * PA, PA degrees (th, ph, mu, la) = {q.degrees()}, {len(q)} terms', flush=True)
ok &= all(nm in ('th', 'ph', 'mu', 'la', '1-la', '1-mu', '1-th-ph') for nm, e in rem)
qd = q.compose(x, (1 - x) * y, m_, l_, ctx=ba.ctx)
ok &= cert(f'PA after the triangle map (H = {H})', qd)
print(f'total {time.time() - t00:.1f} s')
print('PASS' if ok else 'FAIL')
