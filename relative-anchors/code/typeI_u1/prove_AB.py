"""Exact Bernstein certificates for inequalities (A) and (B) (normalised S_T = 1, H = 2/5).

(A)  FA(s, s2, sig) > 0 for s, s2 >= 0, s + s2 <= 1, sig in [0, oo]          =>  ||N||_F <= H   (N = Ghat_tt'/sqrt(s s2))
(B)  FB1(s, sig) > 0 and FB2(s, sig) > 0 for s in [0, 1], sig in [0, oo], except FB2 = 0 at (s, sig) = (1, 0)
     =>  lambda_min(Ghat_tt) > H (1 - s)   on 0 <= s <= 1, sig > 0
(D)  the denominators: E3 > 0 on the same triangle x [0, oo] (the other factors of Den, Dd are obviously positive)

Domain maps (all onto):
  triangle:  s = a, s2 = (1 - a) b, (a, b) in [0,1]^2;
  sigma:     chart I  sig = x in [0, 1];  chart II  sig = 1/x, x in [0, 1], polynomial multiplied by x^deg_sig
             (x = 0 is the limit sig -> oo, i.e. the leading sig-coefficient).
"""
import pickle, time, itertools
import sympy as sp
from math import lcm, comb
from fractions import Fraction as Fr
from bern import from_dict, certify, bernstein_on_box
from polys import s, s2, sg

t00 = time.time()
d = pickle.load(open('polys.pkl', 'rb'))
H = d['H']
print(f"H = {H}", flush=True)


def int_dict(expr, gens):
    P = sp.Poly(expr, *gens)
    den = 1
    for c in P.coeffs():
        den = lcm(den, int(sp.denom(c)))
    return {m: int(c * den) for m, c in P.terms()}, tuple(P.degree_list())


def chart(dct, degs, axis, which):
    """which = 'I': unchanged; 'II': x^deg * p(.., 1/x, ..) (reverse the axis)."""
    if which == 'I':
        return from_dict(dct, degs)
    D = degs[axis]
    nd = {}
    for m, c in dct.items():
        mm = list(m)
        mm[axis] = D - mm[axis]
        nd[tuple(mm)] = c
    return from_dict(nd, degs)


def duffy(dct):
    nd = {}
    for (i, j, k), c in dct.items():
        for r in range(j + 1):
            key = (i + r, j, k)
            nd[key] = nd.get(key, 0) + c * comb(j, r) * (-1) ** r
    nd = {m: c for m, c in nd.items() if c != 0}
    return nd, tuple(max(m[q] for m in nd) for q in range(3))


def root_stats(C):
    """min over K of the true Bernstein coefficient b_K, relative to max |b_K| (root box)."""
    B = bernstein_on_box(C, ((0, 0),) * C.ndim)
    degs = [n - 1 for n in C.shape]
    vals = []
    for K in itertools.product(*[range(n) for n in C.shape]):
        den = 1
        for dj, kj in zip(degs, K):
            den *= comb(dj, kj)
        vals.append(Fr(B[K], den))
    pos = [v for v in vals if v != 0]
    return len(vals), float(min(pos) / max(abs(v) for v in vals)), sum(1 for v in vals if v == 0)


ok_all = True
results = []


def run(name, C, corner=None):
    global ok_all
    t0 = time.time()
    ok, st = certify(C, corner=corner)
    ncoef, rel, nzero = root_stats(C)
    print(f"{name}: {'certified' if ok else 'FAILED'} {st}; {ncoef} Bernstein coefficients on the root box, "
          f"{nzero} zero, smallest positive / largest = {rel:.2e} ({time.time() - t0:.1f} s)", flush=True)
    ok_all &= ok
    results.append((name, ok, st))


# ---------------- (B) ----------------
for name in ('FB1', 'FB2'):
    dct, degs = int_dict(d[name], (s, sg))
    print(f"(B) {name}: degrees (s, sig) = {degs}, {len(dct)} terms", flush=True)
    for which in ('I', 'II'):
        corner = (1, 0) if (name == 'FB2' and which == 'I') else None
        run(f"(B) {name} chart {which}" + (" [corner rule at (s, sig) = (1, 0)]" if corner else ""),
            chart(dct, degs, 1, which), corner)

e = sp.symbols('e')
F2 = sp.Poly(sp.expand(d['FB2'].subs(s, 1 - e)), e, sg)
lin = {m: c for m, c in F2.terms() if sum(m) <= 1}
print(f"(B) FB2 at (s, sig) = (1 - e, sig): constant term {lin.get((0, 0), 0)}, "
      f"linear part {lin.get((1, 0), 0)} e + {lin.get((0, 1), 0)} sig", flush=True)
ok_all &= lin.get((0, 0), 0) == 0 and lin.get((1, 0), 0) > 0 and lin.get((0, 1), 0) > 0

# ---------------- (A) ----------------
dctA, _ = int_dict(d['FA'], (s, s2, sg))
nd, degs2 = duffy(dctA)
print(f"(A) FA after the triangle map: degrees (a, b, x) = {degs2}, {len(nd)} terms", flush=True)
for which in ('I', 'II'):
    run(f"(A) FA chart {which}", chart(nd, degs2, 2, which))

# ---------------- (D) denominators ----------------
dctE, _ = int_dict(d['E3'], (s, s2, sg))
ndE, degsE = duffy(dctE)
for which in ('I', 'II'):
    run(f"(D) E3 chart {which}", chart(ndE, degsE, 2, which))

pickle.dump(results, open('prove_AB_results.pkl', 'wb'))
print(f"total {time.time() - t00:.1f} s")
print("PASS" if ok_all else "FAIL")
