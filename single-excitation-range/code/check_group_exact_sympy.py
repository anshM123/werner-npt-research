"""EXACT certificate for the bath-reduced relative-anchor conditions (README, Lemma 13; (R2) always holds, see the remark after Lemma 13) of a group vector
    p_G = sum_t w_t |1_t>   on k copies, the weights w_t given as rationals or as parameters (rational functions).
Conditions, with A = (M^{(x)k})^Gamma and B = R_{p_G}^Gamma:
  (D)     D_{p_G} > 0                                        (always: D_p >= 2^{-k} |p|^2 I, the rank-one bound);
  (R_inf) Y = A - B compressed to G1_j + G0_{j+1} is PSD,    j = -1..k;
  (R2)    3A + B compressed to G0_j is PSD,                  j = 0..k.
Usage:  python check_group_exact_sympy.py "1,v,u"            (parameters v, u > 0; field QQ(v, u))
        python check_group_exact_sympy.py "1,1,u"            (one parameter u > 0)
        python check_group_exact_sympy.py "1,2/3,5/7" point   (a rational point; also dumps float matrices to point_dump.json)
        python check_group_exact_sympy.py "1,1,u" bscale=1001/1000   (sanity: B replaced by bscale * B in (R_inf); must FAIL)
        python check_group_exact_sympy.py "1,1,u" r2coef=c           (sanity: 3A + B replaced by cA + B in (R2))
Certificate: exact block-diagonality for the torus weights (n_1; n_2 at each copy); in each block symmetric Gaussian
elimination over the field with diagonal pivoting.  A pivot is accepted only if it is certified > 0 on the open
parameter orthant: numerator and denominator each with all coefficients of one sign (and equal signs); in the
one-parameter case alternatively no real root in (0, oo) (exact sympy root isolation) and positive at 1.  An identically
zero remaining diagonal entry must come with an identically zero row.  Then each block is congruent, at every point of
the open orthant, to diag(pivots, 0), hence PSD.  Prints PASS/FAIL."""
import itertools
import json
import sys
import time
import sympy as sp
from sympy import QQ
from sympy.polys.fields import field

t0 = time.time()
spec = sys.argv[1]
opts = sys.argv[2:]
toks = [s.strip() for s in spec.split(',')]
k = len(toks)
names = sorted({s for s in toks if s[0].isalpha()})
point = 'point' in opts
bscale = QQ(1)
r2coef = QQ(3)
for o in opts:
    if o.startswith('r2coef='):
        r2coef = QQ(*map(int, o.split('=')[1].split('/'))) if '/' in o else QQ(int(o.split('=')[1]))
    if o.startswith('bscale='):
        bscale = QQ(*map(int, o.split('=')[1].split('/'))) if '/' in o else QQ(int(o.split('=')[1]))
gens = names if names else ['zdummy']
K, *G = field(','.join(gens), QQ)
G = dict(zip(gens, G))


def tok2(s):
    if s[0].isalpha():
        return G[s]
    if '/' in s:
        a, b = s.split('/')
        return K(QQ(int(a), int(b)))
    return K(QQ(int(s)))


W = [tok2(s) for s in toks]
nvar = len(names)
usym = sp.Symbol(names[0]) if nvar == 1 else None
half = QQ(1, 2)
ZERO, ONE = K(0), K(1)


def m1(a, b, ap, bp):
    return QQ(int(a == ap and b == bp)) - half * int(a == b and ap == bp)


def Mel(row, col, rowp, colp):
    v = QQ(1)
    for s in range(k):
        v *= m1(row[s], col[s], rowp[s], colp[s])
        if v == 0:
            return QQ(0)
    return v


def q(st):
    return sum(x != 0 for x in st)


allst = list(itertools.product(range(3), repeat=k))
e = [tuple(1 if r == t else 0 for r in range(k)) for t in range(k)]
p = {e[t]: W[t] for t in range(k) if W[t] != 0}


def wcol(y):
    return (sum(v == 1 for v in y),) + tuple(int(v == 2) for v in y)


def wstate(a, y):
    return (sum(v == 1 for v in a) + sum(v == 1 for v in y),) + tuple(int(a[s] == 2) + int(y[s] == 2) for s in range(k))


def onesign(poly):
    cs = poly.coeffs()
    return all(c > 0 for c in cs) or all(c < 0 for c in cs)


def no_positive_root(poly):
    P = sp.Poly(poly.as_expr(), usym)
    if P.degree() <= 0:
        return P.as_expr() != 0
    return len([r for r in sp.real_roots(P) if r > 0]) == 0


def val1(poly):
    return sp.Poly(poly.as_expr(), usym).eval(1)


def positive_on_orthant(f):
    if f == 0:
        return False
    num, den = f.numer, f.denom
    if onesign(num) and onesign(den):
        return (num.coeffs()[0] > 0) == (den.coeffs()[0] > 0)
    if nvar == 1 and no_positive_root(num) and no_positive_root(den):
        return (val1(num) > 0) == (val1(den) > 0)
    return False


def inv(Mx):
    n = len(Mx)
    A = [row[:] + [ONE if i == j else ZERO for j in range(n)] for i, row in enumerate(Mx)]
    for c in range(n):
        r = next(i for i in range(c, n) if A[i][c] != 0)
        A[c], A[r] = A[r], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return [row[n:] for row in A]


def leading_pivots_positive(Bm):
    M = [row[:] for row in Bm]
    n = len(M)
    for i in range(n):
        if not positive_on_orthant(M[i][i]):
            return False
        for j in range(i + 1, n):
            if M[j][i] != 0:
                f = M[j][i] / M[i][i]
                for l in range(i, n):
                    if M[i][l] != 0:
                        M[j][l] -= f * M[i][l]
    return True


def psd_certificate(Mx):
    M = [row[:] for row in Mx]
    act = list(range(len(M)))
    npiv = 0
    while act:
        for i in act:
            if M[i][i] == 0 and any(M[i][j] != 0 for j in act):
                return False, None
        cand = [i for i in act if M[i][i] != 0]
        if not cand:
            return True, npiv
        cand.sort(key=lambda i: len(M[i][i].numer.terms()) + len(M[i][i].denom.terms()))
        i = next((c_ for c_ in cand if positive_on_orthant(M[c_][c_])), None)
        if i is None:
            return False, None
        npiv += 1
        rest = [j for j in act if j != i]
        for j in rest:
            if M[j][i] == 0:
                continue
            f = M[j][i] / M[i][i]
            for l in rest:
                if M[i][l] != 0:
                    M[j][l] -= f * M[i][l]
        act = rest
    return True, npiv


# ---- D, inverse by blocks
D = {}
for y in allst:
    for yp in allst:
        if wcol(y) != wcol(yp):
            continue
        v = ZERO
        for a, pa in p.items():
            for b, pb in p.items():
                mv = Mel(a, y, b, yp)
                if mv:
                    v += pa * pb * K(mv)
        D[(y, yp)] = v
# block diagonality of D (exact): entries between different column weights vanish
bdD = all(sum((pa * pb * K(Mel(a, y, b, yp)) for a, pa in p.items() for b, pb in p.items()), ZERO) == 0
          for y in allst for yp in allst if wcol(y) != wcol(yp) and q(y) == q(yp))
blocks = {}
for y in allst:
    blocks.setdefault(wcol(y), []).append(y)
D_ok = bdD
Dinv = {}
for key, ys in blocks.items():
    Bm = [[D[(y, yp)] for yp in ys] for y in ys]
    D_ok &= leading_pivots_positive(Bm)
    Bi = inv(Bm)
    for i, y in enumerate(ys):
        for j, yp in enumerate(ys):
            Dinv[(y, yp)] = Bi[i][j]
print(f'k = {k}, weights ({spec}), parameters {names}: (D) D_pG block diagonal and > 0 (exact): {D_ok}  '
      f'[{time.time() - t0:.1f}s]', flush=True)
rows_needed = [st for st in allst if q(st) <= 1]
MVrow = {}
for b in rows_needed:
    for y in allst:
        for z in allst:
            v = ZERO
            for a, pa in p.items():
                mv = Mel(b, y, a, z)
                if mv:
                    v += K(mv) * pa
            if v != 0:
                MVrow.setdefault((b, y), {})[z] = v
Trow = {}
for key, zs in MVrow.items():
    acc = {}
    for z, v in zs.items():
        for zp in blocks[wcol(z)]:
            dv = Dinv[(z, zp)]
            if dv != 0:
                acc[zp] = acc.get(zp, ZERO) + v * dv
    Trow[key] = acc
Rcache = {}


def Rel_(b, yp, a, y):
    key = (b, yp, a, y)
    if key not in Rcache:
        r1 = Trow.get((b, yp), {})
        r2 = MVrow.get((a, y), {})
        s = ZERO
        for z, v in r1.items():
            v2 = r2.get(z)
            if v2 is not None:
                s += v * v2
        Rcache[key] = s
    return Rcache[key]


def Aent(a, yp, b, y):      # A[(a,y'),(b,y)] = M[(b,y'),(a,y)]
    return K(Mel(b, yp, a, y))


def Bent(a, yp, b, y):      # B[(a,y'),(b,y)] = R[(b,y'),(a,y)]
    return Rel_(b, yp, a, y)


vac = (0,) * k
H = {j: [st for st in allst if q(st) == j] for j in range(k + 1)}
ok = D_ok
dump = {}
for j in range(-1, k + 1):
    S = ([(a, y) for a in H[1] for y in H[j]] if j >= 0 else []) + ([(vac, y) for y in H[j + 1]] if j + 1 <= k else [])
    G0 = [(vac, y) for y in H[j]] if j >= 0 else []
    for name, states, entry in [('R_inf', S, lambda a, yp, b, y: Aent(a, yp, b, y) - bscale * Bent(a, yp, b, y)),
                                ('R2', G0, lambda a, yp, b, y: r2coef * Aent(a, yp, b, y) + Bent(a, yp, b, y))]:
        if not states:
            continue
        tb = {}
        for st in states:
            tb.setdefault(wstate(*st), []).append(st)
        bd = all(entry(a, yp, b, y) == 0 for (a, yp) in states for (b, y) in states if wstate(a, yp) != wstate(b, y))
        allok = bd
        maxsz = rk = 0
        for key, sts in tb.items():
            Mx = [[entry(a, yp, b, y) for (b, y) in sts] for (a, yp) in sts]
            sym = all(Mx[i][l] == Mx[l][i] for i in range(len(sts)) for l in range(len(sts)))
            good, r = psd_certificate(Mx)
            if not (sym and good):
                print('   not certified: block', key, 'size', len(sts), 'symmetric', sym)
            allok &= sym and good
            maxsz = max(maxsz, len(sts))
            rk += r if r else 0
        if point:
            dump[f'{name}_{j}'] = {'states': [[list(a), list(y)] for (a, y) in states],
                                   'M': [[float(sp.Rational(str(entry(a, yp, b, y).numer.coeffs()[0] if entry(a, yp, b, y) != 0 else 0)) /
                                                sp.Rational(str(entry(a, yp, b, y).denom.coeffs()[0])))
                                          if entry(a, yp, b, y) != 0 else 0.0
                                          for (b, y) in states] for (a, yp) in states]}
        ok &= allok
        print(f'  j = {j}: ({name}) {len(states)} states, {len(tb)} torus blocks (largest {maxsz}), block-diagonal {bd}, '
              f'certified PSD: {allok}, rank {rk}  [{time.time() - t0:.1f}s]', flush=True)
if point:
    json.dump({'w': toks, 'mats': dump}, open('point_dump.json', 'w'))
print('PASS' if ok else 'FAIL')
