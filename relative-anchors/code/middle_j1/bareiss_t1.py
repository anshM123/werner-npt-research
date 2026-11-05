"""Leading-minor certificate for the gauge-fixed type I block T = {t} of (R_inf)_1 at a given k (default 5), by
fraction-free (Bareiss) elimination with the denominator factors pulled out.  Adapted from
relative-anchors/code/middle_k4/bareiss_cert.py (same method; NOTES.md there, Section 3).

The block is built exactly over Q(s_0, ..., s_{k-1}) with s_1 = 1 (t = 0; the rescaled entries are homogeneous of
degree 1) by lazyb.LazyBuilder (same entries as ylib.Builder).  Gauge: delete the state (1@1, 2@t) (the kernel vector
p (x) |2@t> is the all-ones vector on the states (1@r, 2@t) and is the only known kernel vector).  Natural order
(row 2@t states, then row 1@r states, then vacuum states).  If the common denominator factors and every leading-minor
core have only positive coefficients, the gauge-fixed block is positive definite on the open orthant (Sylvester), so
the block is >= 0 with kernel exactly p (x) |2@t> for all positive weights.
Usage: python bareiss_t1.py [k]"""
import sys, time, random
from math import lcm
from fractions import Fraction as Fr
import flint
from flintfield import Field
from ylib import Ring
from lazyb import LazyBuilder
from j1lib import block_states, to_ylib_state

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
t0 = time.time()
names = ['s0'] + [f's{i}' for i in range(2, k)]
K = Field(names)
g = K.gens
s = [g[0], K(1)] + list(g[1:])
Bd = LazyBuilder(Ring(s, K(0), K(1)))
st = block_states(k, (0,))
keep = [i for i, (a, y) in enumerate(st) if not (a is not None and a[0] == 1 and a[1] == 1)]
ys = [to_ylib_state(st[i], k) for i in keep]
n = len(keep)
E = [[None] * n for _ in range(n)]
for a in range(n):
    for b in range(a, n):
        E[a][b] = E[b][a] = Bd.Y(ys[a], ys[b])
print(f'k = {k}, T = (0,): gauge-fixed size {n}, built [{time.time() - t0:.1f}s]', flush=True)
ctxZ = flint.fmpz_mpoly_ctx.get(tuple(names), 'lex')
L = 1
for row in E:
    for v in row:
        for c in list(v.n.coeffs()) + list(v.d.coeffs()):
            L = lcm(L, int(flint.fmpq(c).q))


def toZ(p):
    return ctxZ.from_dict({m: int(flint.fmpq(c) * L) for m, c in p.to_dict().items()})


NUMS = [[toZ(E[a][b].n) for b in range(n)] for a in range(n)]
DENS = [[toZ(E[a][b].d) for b in range(n)] for a in range(n)]
den = DENS[0][0]
for a in range(n):
    for b in range(n):
        gg = den.gcd(DENS[a][b])
        den = den * divmod(DENS[a][b], gg)[0]
Nt = [[NUMS[a][b] * divmod(den, DENS[a][b])[0] for b in range(n)] for a in range(n)]
if den.leading_coefficient() < 0:
    den = -den
    Nt = [[-x for x in row] for row in Nt]
cst, fac = den.factor()
FACS = [f for f, e in fac]
nf = len(FACS)
print('  den = ' + str(cst) + ' * ' + ' * '.join(f'[deg {f.total_degree()}, {len(f)} t]^{e}' for f, e in fac), flush=True)
denpos = cst > 0 and all(all(c > 0 for c in f.coeffs()) for f in FACS)
print('  den constant > 0 and den factors with positive coefficients:', denpos, flush=True)


def exact_div(p, q):
    """p / q if q divides p (division with remainder; the plain exact division crashed once on large inputs)."""
    a, r = divmod(p, q)
    return a if r == 0 else None


def strip(p):
    ex = [0] * nf
    if p == 0:
        return p, ex
    for i, f in enumerate(FACS):
        while True:
            q = exact_div(p, f)
            if q is None:
                break
            p = q
            ex[i] += 1
    return p, ex


ZERO = ctxZ.from_dict({})
b = [[None] * n for _ in range(n)]
for a in range(n):
    for c in range(a, n):
        b[a][c] = strip(Nt[a][c])
leading = [b[0][0]]
prev = (ctxZ.from_dict({(0,) * len(names): 1}), [0] * nf)
for p in range(n - 1):
    cp, ep = b[p][p]
    if cp == 0:
        print('  zero pivot', p, flush=True)
        break
    cprev, eprev = prev
    for i in range(p + 1, n):
        ci, ei = b[p][i]
        for l in range(i, n):
            cl, el = b[p][l]
            cil, eil = b[i][l]
            e1 = [x + y for x, y in zip(ep, eil)]
            e2 = [x + y for x, y in zip(ei, el)]
            if cil == 0:
                gq = e2
            elif ci == 0 or cl == 0:
                gq = e1
            else:
                gq = [min(x, y) for x, y in zip(e1, e2)]
            t = ZERO
            if cil != 0:
                m1 = cp * cil
                for f, x, g2 in zip(FACS, e1, gq):
                    if x > g2:
                        m1 = m1 * f ** (x - g2)
                t = t + m1
            if ci != 0 and cl != 0:
                m2 = ci * cl
                for f, x, g2 in zip(FACS, e2, gq):
                    if x > g2:
                        m2 = m2 * f ** (x - g2)
                t = t - m2
            if t == 0:
                b[i][l] = (ZERO, [0] * nf)
                continue
            q = exact_div(t, cprev)
            assert q is not None, 'inexact division by previous core'
            ex = [g2 - x for g2, x in zip(gq, eprev)]
            for kk, f in enumerate(FACS):
                while ex[kk] < 0:
                    q2 = exact_div(q, f)
                    assert q2 is not None, 'inexact division by a factor'
                    q = q2
                    ex[kk] += 1
            c2, e3 = strip(q)
            b[i][l] = (c2, [x + y for x, y in zip(ex, e3)])
    prev = b[p][p]
    leading.append(b[p + 1][p + 1])
    lc = b[p + 1][p + 1][0]
    print(f'  step {p + 1}: D_{p + 2} core deg {lc.total_degree()}, {len(lc)} terms, negative coeffs '
          f'{sum(1 for c in lc.coeffs() if c < 0)} [{time.time() - t0:.1f}s]', flush=True)
    for l in range(p, n):
        b[p][l] = None if l != p else b[p][l]
random.seed(11)
pt = [random.randint(2, 30) for _ in names]
Nv = [[Fr(int(Nt[a][c](*pt))) for c in range(n)] for a in range(n)]


def detF(A):
    A = [row[:] for row in A]
    m = len(A)
    d = Fr(1)
    for p in range(m):
        piv = next((i for i in range(p, m) if A[i][p] != 0), None)
        if piv is None:
            return Fr(0)
        if piv != p:
            A[p], A[piv] = A[piv], A[p]
            d = -d
        d *= A[p][p]
        for i in range(p + 1, m):
            f = A[i][p] / A[p][p]
            for l in range(p, m):
                A[i][l] -= f * A[p][l]
    return d


ok = denpos
allpos = True
for r, (c, ex) in enumerate(leading, 1):
    val = Fr(int(c(*pt)))
    for f, e in zip(FACS, ex):
        val *= Fr(int(f(*pt))) ** e
    good = val == detF([row[:r] for row in Nv[:r]])
    pos = all(x > 0 for x in c.coeffs())
    ok &= good
    allpos &= pos
    print(f'  leading minor {r}: point check {good}; core deg {c.total_degree()}, {len(c)} terms, '
          f'{"all coefficients positive" if pos else "MIXED"}', flush=True)
ok &= allpos and len(leading) == n
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
