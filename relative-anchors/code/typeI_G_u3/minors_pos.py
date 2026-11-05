"""Proof step: G > 0 for one pair (u, m) from the closed form closed_u{u}_m{m}.pkl (closed_gen.py), exact.

G_z (z-coordinates, order rho_0, z_{0,h}..., rho_1, ...) is congruent to G by an invertible diagonal matrix, so G > 0
iff G_z > 0.  Write G_z = Nt/den with one common denominator den = cst * prod_i f_i^{E_i} (f_i irreducible, integer
coefficients after scaling).  For r = 1..n the minors M(rows 0..r-1 | J) (J a set of r columns) are computed by
expansion along row r - 1 from those of level r - 1, each stored as prod_i f_i^{x_i} * core with the factors f_i pulled
out by exact division.  The leading minors of G_z are
        minor_r(G_z) = core_r * prod_i f_i^{x_ir - r E_i} / cst^r.
If cst > 0 and every f_i and every core_r has only positive coefficients, every leading minor is positive on the open
orthant (all s_i > 0), so G_z > 0 there by Sylvester's criterion.  The script checks exactly: the expansions (every
leading minor against a direct determinant of Nt at a rational point), cst > 0, and the coefficient signs.  It prints
PASS only if all of these hold.  Usage: python minors_pos.py u m [dehom]
With 'dehom' the closed form computed with s0 = 1 is used (closed_gen.py ... dehom); positivity of all coefficients
of the dehomogenised polynomials is the same as for their homogenisations, and G_z(lambda s) is congruent to G_z(s)
by a positive diagonal matrix, so the conclusion is the same."""
import pickle, sys, time, flint, random, itertools
from math import lcm
from fractions import Fraction as Fr

u, m = int(sys.argv[1]), int(sys.argv[2])
dehom = len(sys.argv) > 3 and sys.argv[3] == 'dehom'
fn = f'closed_u{u}_m{m}{"_dehom" if dehom else ""}.pkl'
t0 = time.time()
D = pickle.load(open(fn, 'rb'))
names = tuple(D['names'])
ctxZ = flint.fmpz_mpoly_ctx.get(names, 'lex')
G = D['G']
n = len(G)
L = 1
for row in G:
    for (n_, d_) in row:
        for c in list(n_.values()) + list(d_.values()):
            L = lcm(L, int(c.q))


def toZ(dct):
    return ctxZ.from_dict({mm: int(c * L) for mm, c in dct.items()})


NUMS = [[toZ(G[i][j][0]) for j in range(n)] for i in range(n)]
DENS = [[toZ(G[i][j][1]) for j in range(n)] for i in range(n)]
den = DENS[0][0]
for i in range(n):
    for j in range(n):
        g = den.gcd(DENS[i][j])
        den = den * (DENS[i][j] / g) if not (DENS[i][j] / g).is_one() else den
for i in range(n):
    for j in range(n):
        q, r = divmod(den, DENS[i][j])
        assert r == 0
Nt = [[NUMS[i][j] * divmod(den, DENS[i][j])[0] for j in range(n)] for i in range(n)]
if den.leading_coefficient() < 0:
    den = -den
    Nt = [[-x for x in row] for row in Nt]
assert all(Nt[i][j] == Nt[j][i] for i in range(n) for j in range(n))
cst, fac = den.factor()
FACS = [f for f, e in fac]
E = [e for f, e in fac]
print(f'u = {u}, m = {m}, n = {n}: den = {cst} * ' + ' * '.join(f'[deg {f.total_degree()}, {len(f)} t]^{e}' for f, e in fac), flush=True)
print('  den factors with all coefficients positive:', [all(c > 0 for c in f.to_dict().values()) for f in FACS], flush=True)


def divout(p, q):
    e = 0
    while True:
        qq, r = divmod(p, q)
        if r != 0:
            return p, e
        p = qq
        e += 1


def strip(p):
    ex = []
    for f in FACS:
        p, e = divout(p, f)
        ex.append(e)
    return p, ex


level = {(): (ctxZ.from_dict({(0,) * len(names): 1}), [0] * len(FACS))}
leading = []
for r in range(1, n + 1):
    new = {}
    need_cols = range(n)
    for J in itertools.combinations(range(n), r):
        if J[-1] > max(r - 1, 0) + (n - r):   # all subsets are needed for later levels
            pass
        items = []
        for pos, j in enumerate(J):
            sub = tuple(x for x in J if x != j)
            items.append(((-1) ** (r - 1 + pos), j, level[sub]))
        g = [min(it[2][1][i] for it in items) for i in range(len(FACS))]
        tot = 0 * Nt[0][0]
        for sgn, j, (c, ex) in items:
            p = c
            for f, e, gi in zip(FACS, ex, g):
                if e > gi:
                    p = p * f ** (e - gi)
            tot = tot + sgn * Nt[r - 1][j] * p
        c2, e2 = strip(tot)
        new[J] = (c2, [a + b for a, b in zip(g, e2)])
    level = new
    lead = level[tuple(range(r))]
    leading.append(lead)
    print(f'  level {r} done: {len(level)} minors ({time.time() - t0:.1f} s); leading core degree {lead[0].total_degree()}, '
          f'{len(lead[0])} terms, exponents {lead[1]}', flush=True)
    if r < n:
        level = {J: v for J, v in level.items()}
# check at a rational point
random.seed(5)
pt = [random.randint(1, 25) for _ in names]
Nv = [[Fr(int(Nt[i][j](*pt))) for j in range(n)] for i in range(n)]


def detF(A):
    A = [row[:] for row in A]
    k = len(A)
    det_ = Fr(1)
    for p in range(k):
        piv = next((i for i in range(p, k) if A[i][p] != 0), None)
        if piv is None:
            return Fr(0)
        if piv != p:
            A[p], A[piv] = A[piv], A[p]
            det_ = -det_
        det_ *= A[p][p]
        for i in range(p + 1, k):
            f = A[i][p] / A[p][p]
            for j in range(p, k):
                A[i][j] -= f * A[p][j]
    return det_


ok = True
allpos = True
for r, (c, ex) in enumerate(leading, 1):
    val = Fr(int(c(*pt)))
    for f, e in zip(FACS, ex):
        val *= Fr(int(f(*pt))) ** e
    good = val == detF([row[:r] for row in Nv[:r]])
    ok &= good
    coeffs = list(c.to_dict().values())
    pos = all(x > 0 for x in coeffs)
    neg = all(x < 0 for x in coeffs)
    allpos &= pos
    print(f'  leading minor {r}: check {good}; core coefficients: {"all positive" if pos else ("all negative" if neg else "MIXED")} '
          f'({sum(1 for x in coeffs if x < 0)} negative of {len(coeffs)})', flush=True)
pickle.dump(dict(names=names, den_factors=[(f.to_dict(), e) for f, e in fac], den_const=int(cst),
                 cores=[(c.to_dict(), ex) for c, ex in leading]), open(f'minors_u{u}_m{m}{"_dehom" if dehom else ""}.pkl', 'wb'))
denpos = cst > 0 and all(all(c > 0 for c in f.to_dict().values()) for f in FACS)
print(f'expansion checks: {ok}; den constant > 0 and all den factors with positive coefficients: {denpos}; '
      f'all leading-minor cores with positive coefficients: {allpos}')
if ok and denpos and allpos:
    print(f'CONCLUSION: every leading principal minor of G_z is positive on the open orthant, so G > 0 for u = {u}, m = {m}, all positive weights')
ok = ok and denpos and allpos
print(f'time {time.time() - t0:.1f} s')
print('PASS' if ok else 'FAIL')
