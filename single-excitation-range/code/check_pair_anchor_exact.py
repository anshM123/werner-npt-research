"""EXACT verification (sympy, rational arithmetic in one parameter u) of the two-copy facts behind the pair theorem.

p_u = |10> + u|01> (unnormalised; T_p and the sign of D depend only on the ray), u in [0, oo).  Every two-copy vector of
the form x (x) |0> + |0> (x) y with x, y orthogonal to |0> is, after local unitaries fixing |0> and a global phase,
a multiple of some p_u or equal to |01>.  We prove, for every u >= 0:
  (F1) D_u = V^* M^{(x)2} V > 0;
  (F2) A - R^Gamma >= 0          (T <= I);
  (F3) R^Gamma + 3A >= 0          (T >= -3, so that appending product copies keeps T <= I).
Method: A and R^Gamma are block diagonal for the torus weights (n_1, n_{2,copy 0}, n_{2,copy 1}) (checked exactly).
On each block H, H >= 0 iff every elementary symmetric function E_k(H) (sum of the k x k principal minors) is >= 0.
After multiplying by a positive denominator each E_k is a polynomial in u with rational coefficients; its
nonnegativity on [0, oo) is certified by exact real-root isolation plus exact evaluation at a rational point in every
gap between consecutive roots (and beyond the last one)."""
import sympy as sp
import itertools
import time

u, x = sp.symbols('u x')
d = 3
t0 = time.time()


def M1(a, b, ap, bp):
    return sp.Integer(int(a == ap and b == bp)) - sp.Rational(1, 2) * int(a == b and ap == bp)


rows = list(itertools.product(range(d), repeat=2))          # (a0, a1)
cols = rows
full = [(r, c) for r in rows for c in cols]                  # vec index ((a0,a1),(b0,b1))
fidx = {st: i for i, st in enumerate(full)}


def M2(st, stp):
    (a0, a1), (b0, b1) = st
    (ap0, ap1), (bp0, bp1) = stp
    return M1(a0, b0, ap0, bp0) * M1(a1, b1, ap1, bp1)


p = {(1, 0): sp.Integer(1), (0, 1): u}
MV = sp.zeros(81, 9)
for i, st in enumerate(full):
    for j, y in enumerate(cols):
        MV[i, j] = sum(M2(st, (ap, y)) * pv for ap, pv in p.items())
D = sp.zeros(9, 9)
for i, y in enumerate(cols):
    for j, yp in enumerate(cols):
        D[i, j] = sp.expand(sum(pa * pb * M2((a, y), (b, yp)) for a, pa in p.items() for b, pb in p.items()))


def wcol(y):
    return (sum(v == 1 for v in y), int(y[0] == 2), int(y[1] == 2))


def nonneg_on_halfline(f, strict=False):
    """exact: f(u) >= 0 (> 0 if strict) for all u >= 0; f a polynomial in u with rational coefficients."""
    P = sp.Poly(sp.expand(f), u, domain='QQ')
    if P.is_zero:
        return not strict
    ints = P.intervals(inf=0, eps=sp.Rational(1, 10 ** 6))
    ints = sorted((sp.Rational(a), sp.Rational(b)) for (a, b), m in ints)
    if strict and ints:
        return False
    pts = [sp.Integer(0)]
    for k in range(len(ints)):
        lo = ints[k][1]
        if k + 1 < len(ints):
            pts.append((lo + ints[k + 1][0]) / 2)
        else:
            pts.append(lo + 1)
    if not ints:
        pts.append(sp.Integer(1))
    if strict:
        return all(P.eval(q) > 0 for q in pts)
    return all(P.eval(q) >= 0 for q in pts)


ok = True
dg = {}
for j, y in enumerate(cols):
    dg.setdefault(wcol(y), []).append(j)
bd = all(D[i, j] == 0 for i in range(9) for j in range(9) if wcol(cols[i]) != wcol(cols[j]))
Dok = bd
Dinfo = []
for k, js in dg.items():
    B = D.extract(js, js)
    Dinfo.append((k, [sp.factor(B.extract(list(range(r)), list(range(r))).det()) for r in range(1, len(js) + 1)]))
    for r in range(1, len(js) + 1):
        if not nonneg_on_halfline(B.extract(list(range(r)), list(range(r))).det(), strict=True):
            Dok = False
print('(F1) D_u block diagonal and > 0 for all u >= 0 (exact, leading principal minors):', Dok)
for k, mins in Dinfo:
    print('     weight', k, 'leading minors', mins)
ok &= Dok
Dinv = sp.zeros(9, 9)
for k, js in dg.items():
    Binv = D.extract(js, js).inv().applyfunc(sp.cancel)
    for a_, ja in enumerate(js):
        for b_, jb in enumerate(js):
            Dinv[ja, jb] = Binv[a_, b_]
R = (MV * Dinv * MV.T).applyfunc(sp.cancel)


def PT(W):
    out = sp.zeros(81, 81)
    for i, (a, yp) in enumerate(full):
        for j, (b, y) in enumerate(full):
            out[i, j] = W[fidx[(b, yp)], fidx[(a, y)]]
    return out


Mfull = sp.Matrix(81, 81, lambda i, j: M2(full[i], full[j]))
A = PT(Mfull)
RG = PT(R)


def wt(st):
    (a0, a1), (b0, b1) = st
    return (sum(v == 1 for v in (a0, a1, b0, b1)), int(a0 == 2) + int(b0 == 2), int(a1 == 2) + int(b1 == 2))


groups = {}
for i, st in enumerate(full):
    groups.setdefault(wt(st), []).append(i)
blockdiag = all(A[i, j] == 0 and RG[i, j] == 0 for i in range(81) for j in range(81) if wt(full[i]) != wt(full[j]))
print('A and R^Gamma block diagonal for the torus weights (exact):', blockdiag,
      f'({len(groups)} blocks, largest {max(len(v) for v in groups.values())})  [{time.time() - t0:.0f}s]')
ok &= blockdiag


def psd_on_halfline(H):
    den = sp.Integer(1)
    for e in H:
        den = sp.lcm(den, sp.denom(sp.together(e)))
    if not nonneg_on_halfline(den, strict=True):
        den = -den
    assert nonneg_on_halfline(den, strict=True)
    P = (H * den).applyfunc(lambda e: sp.expand(sp.cancel(e)))
    cp_ = P.charpoly(x).all_coeffs()
    return all(nonneg_on_halfline((-1) ** k * cp_[k]) for k in range(1, len(cp_)))


res2 = True
res3 = True
for k in sorted(groups):
    idx = groups[k]
    Ab = A.extract(idx, idx)
    Rb = RG.extract(idx, idx)
    r2 = psd_on_halfline(Ab - Rb)
    r3 = psd_on_halfline(Rb + 3 * Ab)
    res2 &= r2
    res3 &= r3
    print(f'  block {k} size {len(idx)}: A - R^G >= 0: {r2};  R^G + 3A >= 0: {r3}   [{time.time() - t0:.0f}s]')
print('(F2) A - R^Gamma >= 0 for all u >= 0 (exact):', res2)
print('(F3) R^Gamma + 3A >= 0 for all u >= 0 (exact):', res3)
ok &= res2 and res3
print('PASS' if ok else 'FAIL')
