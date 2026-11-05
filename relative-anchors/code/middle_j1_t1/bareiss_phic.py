"""Leading-minor certificate, for all positive weights at a fixed k, for the inverse-free congruent form
Phi' = diag(A1, A2) Phi diag(A1, A2) of the P-first reduction (Lemma P), by fraction-free (Bareiss) elimination with
the denominator factors pulled out (the method of middle_j1/bareiss_t1.py, applied to Phi' of size 2k - 1 instead of
the gauge-fixed block of size 3k - 2).

python bareiss_phic.py k                       (all weights at this k)
python bareiss_phic.py levels l off fix which [F|C]   (a level family, e.g. levels 3 1,1,2 1,1,0 inv F: two single
                    U-copies and u_3 = 2 + x_3 equal ones; which = inv, or stdN for level N; F: the matrix Phi of
                    Lemma P, C (default): Phi')
Phi' is built by lev.py (for k: l = k - 1 levels of size one; p = 1; rescaled frame).  For a level family the matrix
is the invariant (or standard) part of D Phi' D, and the identity Phi' = diag(A1, A2) Phi diag(A1, A2) is checked in
the operator algebra of lev.py, exactly (Phi itself is checked against j1lib by check_lev.py).  N = den Phi'
with den the lcm of the entry denominators; every irreducible factor of den must have only positive coefficients.
Bareiss elimination over Z[q] gives D_r(N) = core_r * prod f_i^{e_ir} exactly.  If every core has only positive
coefficients, every leading principal minor of Phi' is positive on the open orthant, so Phi' > 0, Phi > 0, and the
block is >= 0 with kernel exactly p (x) |2@t> for every positive weight vector (Lemma P).  Every D_r is compared with a
direct determinant at a random integer point.  Prints PASS if all hold."""
import sys, time, random
from math import lcm
from fractions import Fraction as Fr
import flint
from lev import build

t0 = time.time()
if sys.argv[1] == 'levels':
    l = int(sys.argv[2])
    offs = [int(v) for v in sys.argv[3].split(',')]
    fix = [v == '1' for v in sys.argv[4].split(',')]
    which = sys.argv[5]
    form = sys.argv[6] if len(sys.argv) > 6 else 'C'
    if form == 'F':      # Phi itself (checked against j1lib by check_lev.py)
        S = build(l, offs, 1, want_phi=True, fixed=fix, p_one=True)
        E = S['Finv'] if which == 'inv' else S['Fstd'][int(which[3:]) - 1]
        cong_ok = True
    else:                # Phi', with the congruence identity checked
        S = build(l, offs, 1, want_phi=True, fixed=fix, p_one=True, want_phic=True, check_cong=True)
        E = S['Cinv'] if which == 'inv' else S['Cstd'][int(which[3:]) - 1]
        cong_ok = S['cong_ok']
    label = f'levels {offs}, fixed {fix}, {which}, form {form}'
    print(f'{label}: congruence identity (form C only): {cong_ok}', flush=True)
else:
    k = int(sys.argv[1])
    u = k - 1
    S = build(u, [1] * u, 1, want_phi=False, fixed=[True] * u, p_one=True, want_phic=True)
    E = S['Cinv']
    label = f'k = {k}, form C (the congruent form)'
    cong_ok = True
n = len(E)
names = S['names']
print(f'{label}: matrix of size {n} over Q({",".join(names)}) (p = 1), built [{time.time() - t0:.1f}s]', flush=True)
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
xi = [i for i, nm in enumerate(names) if nm.startswith('x')]


def xfree(f):
    """f has a monomial free of the x-variables (so a polynomial with positive coefficients is > 0 at x >= 0)."""
    return any(all(m[i] == 0 for i in xi) for m in f.monoms())


denpos = cst > 0 and all(all(c > 0 for c in f.coeffs()) and xfree(f) for f in FACS)
print('  den constant > 0 and den factors with positive coefficients:', denpos, flush=True)


def exact_div(p, q):
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
lc = b[0][0][0]
print(f'  D_1 core deg {lc.total_degree()}, {len(lc)} terms, negative coeffs {sum(1 for c in lc.coeffs() if c < 0)}',
      flush=True)
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


ok = denpos and cong_ok and all(E[a][b] == E[b][a] for a in range(n) for b in range(n))
allpos = True
for r, (c, ex) in enumerate(leading, 1):
    val = Fr(int(c(*pt)))
    for f, e in zip(FACS, ex):
        val *= Fr(int(f(*pt))) ** e
    good = val == detF([row[:r] for row in Nv[:r]])
    pos = all(x > 0 for x in c.coeffs()) and xfree(c)
    ok &= good
    allpos &= pos
    print(f'  leading minor {r}: point check {good}; core deg {c.total_degree()}, {len(c)} terms, '
          f'{"all coefficients positive" if pos else "MIXED"}', flush=True)
ok &= allpos and len(leading) == n
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
