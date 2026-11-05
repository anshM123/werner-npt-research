"""Negative controls and consistency checks for the Bernstein certificates of prove_AB.py.
(1) Bernstein arrays of the actual certificates: vertex coefficients equal the polynomial values at the vertices
    (exact), reconstruction at random rational interior points equals direct evaluation (exact); min coefficients.
(2) Controls: with H = 7/25 (< sup ||N||_F = sqrt(26)/18 = 0.28327...) the (A) certificate must FAIL and FA must be
    negative at the point (s, s2, sig) = (1/2, 1/2, 0); with H = 4/5 (> inf lambda_min/(1-s) = 3/4) the (B)
    certificate must FAIL and FB2 must be negative near (s, sig) = (0, 0)."""
import pickle, random, time
from fractions import Fraction as Fr
from math import comb, lcm
import sympy as sp
from bern import from_dict, bernstein_on_box, certify
from polys import build, s, s2, sg

t0 = time.time()
ok = True
def int_dict(expr, gens):
    P = sp.Poly(expr, *gens); den = 1
    for c in P.coeffs(): den = lcm(den, int(sp.denom(c)))
    return {m: int(c * den) for m, c in P.terms()}, P.degree_list()
def duffy(dctA):
    nd = {}
    for (i, j, k), c in dctA.items():
        for r in range(j + 1):
            key = (i + r, j, k); nd[key] = nd.get(key, 0) + c * comb(j, r) * (-1) ** r
    nd = {m: c for m, c in nd.items() if c != 0}
    return nd, tuple(max(m[q] for m in nd) for q in range(3))
def evald(dct, x):
    tot = Fr(0)
    for m, c in dct.items():
        v = Fr(c)
        for xi, e in zip(x, m):
            if e: v *= Fr(xi) ** e
        tot += v
    return tot
d = pickle.load(open('polys.pkl', 'rb'))
dA, _ = int_dict(d['FA'], (s, s2, sg)); nd, degs = duffy(dA)
C = from_dict(nd, degs); B = bernstein_on_box(C, ((0, 0),) * 3)
# (1) vertices
for v in [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(1,1,1),(1,0,1)]:
    K = tuple(degs[j] if v[j] else 0 for j in range(3))
    ok &= B[K] == evald(nd, v)
# reconstruction
random.seed(3)
for _ in range(3):
    x = [Fr(random.randint(1, 30), 31) for _ in range(3)]
    val = Fr(0)
    for K in [(i, j, k) for i in range(degs[0]+1) for j in range(degs[1]+1) for k in range(degs[2]+1)]:
        bb = Fr(B[K])
        for q in range(3):
            bb *= x[q] ** K[q] * (1 - x[q]) ** (degs[q] - K[q])
        val += bb
    ok &= val == evald(nd, x)
# the scaled coefficients beta_K = prod C(d,k) b_K; report the min of the true b_K relative to the max
bmin = min(Fr(B[K], comb(degs[0], K[0]) * comb(degs[1], K[1]) * comb(degs[2], K[2])) for K in __import__('itertools').product(*[range(dd+1) for dd in degs]))
bmax = max(abs(Fr(B[K], comb(degs[0], K[0]) * comb(degs[1], K[1]) * comb(degs[2], K[2]))) for K in __import__('itertools').product(*[range(dd+1) for dd in degs]))
print(f"(1) FA chart I: vertex values and interior reconstruction exact: {ok}; min/max Bernstein coefficient {float(bmin/bmax):.3e}", flush=True)
# also the Duffy map itself: nd(a, b, x) == FA(a, (1-a) b, x)
for _ in range(3):
    a_, b_, x_ = [Fr(random.randint(0, 30), 30) for _ in range(3)]
    ok &= evald(nd, (a_, b_, x_)) == evald(dA, (a_, (1 - a_) * b_, x_))
print("    Duffy map exact on random points:", ok, flush=True)
# (2) controls
for Hc, name in ((sp.Rational(7, 25), 'FA'), (sp.Rational(4, 5), 'FB2')):
    dc = build(H=Hc, verbose=False)
    if name == 'FA':
        dct, _ = int_dict(dc['FA'], (s, s2, sg)); ndc, degc = duffy(dct)
        okc, st = certify(from_dict(ndc, degc), max_boxes=3000)
        neg = evald(dct, (Fr(1, 2), Fr(1, 2), 0))
    else:
        dct, degc = int_dict(dc['FB2'], (s, sg))
        okc, st = certify(from_dict(dct, degc), corner=(1, 0), max_boxes=3000)
        neg = evald(dct, (Fr(1, 100), Fr(1, 100)))
    print(f"(2) control H = {Hc}, {name}: certificate {'FAILED as it must' if not okc else 'PASSED (BAD)'} ({st.get('reason')}); "
          f"value at the test point {float(neg):.4e} (must be < 0)", flush=True)
    ok &= (not okc) and neg < 0
print(f"time {time.time() - t0:.1f} s")
print("PASS" if ok else "FAIL")
