"""Condition (i) of Lemma F for all weights (NOTES.md, Section 12): D(s_h) > 0.

D (D_natural.pkl) depends on (s, sigma, S) and on g = sum s_k/d_k, N = sum s_k/d_k^2, g0 = sum s_k/d0_k,
N0 = sum s_k/d0_k^2, IG = 1/(2 + g), IGT = 1/(2 + sigma/d0_t + g0).
Part A: N enters only D_aa, D_cc and N0 only D_zz, with coefficients -kappa, kappa >= 0 (one-signed coefficients).
        Since N <= g/c and N0 <= g0/c0 (c = S/2 + sigma <= d_k, c0 = (sigma + S)/2 <= d0_k), D >= Dt := D(N = g/c,
        N0 = g0/c0) in the Loewner order.
Part B: Dt > 0 on the region, normalised K = sigma + S/2 = 1: t = sigma in (0, 1), S = 2(1 - t), s = S xi with
        xi in (0, 1], g = s/(1+s) + th (S - s), g0 = s/(s+b) + th0 (S - s)/b, th, th0 in [0, 1], b = 1 - t/2.
        (g lies in this range because g' = sum_{k != h} s_k/(1 + s_k) is in [0, S - s]; likewise g0.)
        The three leading principal minors of Dt are certified positive: every irreducible factor of numerator and
        denominator is shown one-signed on the domain (positive coefficients, or Bernstein coefficients on
        [0, 1]^4 after removing monomial factors).  Output: logs only."""
import pickle, time, sys
import flint
from flintfield import Field
import bern

t0 = time.time()
P = pickle.load(open('D_natural.pkl', 'rb'))
names = P['names']
iN, iN0 = names.index('N'), names.index('N0')
ok = True
# ---------------- Part A
for key, idx in ((('ah', 'ah'), iN), (('ch', 'ch'), iN), (('zh', 'zh'), iN0)):
    num, den = P['D'][key]
    P1 = {m: c for m, c in num.items() if m[idx] == 1}
    others = [m for m in num if m[idx] > 1]
    neg_ok = all(c < 0 for c in P1.values()) and not others
    den_ok = all(c > 0 for c in den.values())
    print(f'Part A {key}: coefficient of {names[idx]} = -(polynomial with {len(P1)} positive coefficients)/(positive '
          f'denominator): {neg_ok and den_ok}', flush=True)
    ok &= neg_ok and den_ok
for key, (num, den) in P['D'].items():
    if key not in (('ah', 'ah'), ('ch', 'ch')):
        ok &= all(m[iN] == 0 for m in num) and all(m[iN] == 0 for m in den)
    if key != ('zh', 'zh'):
        ok &= all(m[iN0] == 0 for m in num) and all(m[iN0] == 0 for m in den)
print('Part A: N only in D_aa, D_cc and N0 only in D_zz:', ok, flush=True)
# ---------------- Part B
Kf = Field(['t', 'xi', 'th', 'th0'])
t, xi, th, th0 = Kf.gens
one = Kf(1)
sg = t
S = 2 * (one - t)
s = S * xi
b = one - t / 2
c = one
c0 = b
g = s / (one + s) + th * (S - s)
g0 = s / (s + b) + th0 * (S - s) / b
IG = one / (2 + g)
IGT = one / (2 + t / (one + t / 2) + g0)
vals = dict(s=s, sigma=sg, S=S, IG=IG, IGT=IGT, g=g, N=g / c, g0=g0, N0=g0 / c0, u=one)
vv = [vals[n] for n in names]


def pw(v, e):
    r = one
    for _ in range(int(e)):
        r = r * v
    return r


def evalpoly(d):
    tot = Kf(0)
    for m, c_ in d.items():
        term = Kf(c_)
        for v, e in zip(vv, m):
            if e:
                term = term * pw(v, e)
        tot = tot + term
    return tot


Dt = {}
for key, (num, den) in P['D'].items():
    Dt[key] = evalpoly(num) / evalpoly(den)
print(f'Part B: Dt built ({time.time() - t0:.1f} s)', flush=True)
A = [[Dt[('ah', 'ah')], Dt[('ah', 'ch')] / 2, Dt[('ah', 'zh')] / 2],
     [Dt[('ah', 'ch')] / 2, Dt[('ch', 'ch')], Dt[('ch', 'zh')] / 2],
     [Dt[('ah', 'zh')] / 2, Dt[('ch', 'zh')] / 2, Dt[('zh', 'zh')]]]
m1 = A[0][0]
m2 = A[0][0] * A[1][1] - A[0][1] * A[1][0]
m3 = (A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
      + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]))
print(f'minors built ({time.time() - t0:.1f} s)', flush=True)


def factor_sign(f):
    """sign of the polynomial f (fmpq_mpoly in t, xi, th, th0) on t in (0,1), xi in (0,1], th, th0 in [0,1]:
    +1 or -1 if certified, 0 if not.  Monomial factors t, xi are allowed; 1 - t is allowed."""
    dct = f.to_dict()
    cs = list(dct.values())
    if all(cc > 0 for cc in cs) and any(m[2] == 0 and m[3] == 0 for m in dct):
        return 1, 'positive coefficients'
    if all(cc < 0 for cc in cs) and any(m[2] == 0 and m[3] == 0 for m in dct):
        return -1, 'negative coefficients'
    # Bernstein on [0,1]^4 (integer coefficients)
    from math import lcm
    L = 1
    for cc in cs:
        L = lcm(L, int(cc.q))
    degs = [max(m[j] for m in dct) for j in range(4)]
    C = bern.zeros(tuple(d + 1 for d in degs))
    for m, cc in dct.items():
        C[tuple(m)] += int(cc.p) * (L // int(cc.q))
    for sgn in (1, -1):
        okb, st = bern.certify(sgn * C, max_boxes=4000)
        if okb:
            return sgn, f'Bernstein on [0,1]^4, {st}'
    # try with the corner rule at the vertices where it may vanish
    for corner in [(1, 1, 0, 0), (1, 1, 1, 1), (1, 0, 0, 0), (0, 0, 0, 0), (1, 1, 1, 0), (1, 1, 0, 1)]:
        for sgn in (1, -1):
            okb, st = bern.certify(sgn * C, max_boxes=4000, corner=corner)
            if okb:
                return sgn, f'Bernstein with corner rule at {corner}, {st}'
    return 0, 'not certified'


for name, m in (('m1', m1), ('m2', m2), ('m3', m3)):
    sgn_total = 1
    for part, poly in (('numerator', m.n), ('denominator', m.d)):
        cst, facs = poly.factor()
        sgn_part = 1 if cst > 0 else -1
        for f, e in facs:
            sg_f, how = factor_sign(f)
            dct = f.to_dict()
            degs = [max(mm[j] for mm in dct) for j in range(4)]
            print(f'  {name} {part}: factor of degrees {degs}, {len(dct)} terms, exponent {e}: sign {sg_f} ({how})',
                  flush=True)
            if sg_f == 0:
                sgn_part = 0
            else:
                sgn_part *= sg_f ** e
        sgn_total *= sgn_part
    print(f'{name}: certified sign {sgn_total} ({time.time() - t0:.1f} s)', flush=True)
    ok &= sgn_total == 1
print('PASS' if ok else 'FAIL')
