"""Certificate for the type I blocks with |T| = 2 of (R_inf)_1 (NOTES.md, Theorem D1): the 8 reduced inequalities.

Homogeneous coordinates: al = x_t, be = x_t', e = sum of the other x_u (al + be + e = |s| = 1 after scaling; all
quantities below are homogeneous of degree 0 in (al, be, e)).  lam in [0, 1] parametrises Z = F3/e between
Z0 = 1/(al + be + 2e) and Z1 = 1/(al + be).  For each choice (lo/hi) of F0, F1, F2:
    F_i^hi = e Z/(1 + d_i Z),   F_i^lo = e[(1 - lam) h_i(Z0) + lam h_i(Z1)],   h_i(z) = z/(1 + d_i z),
    d_0 = al + be, d_1 = al, d_2 = be,
and g = 4 - c^2 - G3 (A^2 + B^2 - c A B) with
    c = kap/(2(G0 + F0/2)), A = a'/(K1 + F1/2), B = b'/(K2 + F2/2), G3 = 1 + e Z/2.
Claim checked here (exact, python-flint): g = N/D with D a product of factors with only positive coefficients (in
al, be, e and in the Bernstein basis of lam), and N has, in the Bernstein basis of lam of its degree, coefficients that
are polynomials in (al, be, e) with only nonnegative coefficients, the first and the last of them nonzero.  Hence
g > 0 for al, be, e > 0 and every lam in [0, 1].  Also checked: c < 2 in the same sense (2 - c has positive numerator
and denominator).
Usage: python cert_t2.py"""
import itertools, time
import flint
from flintfield import Field

t0 = time.time()
K = Field(['al', 'be', 'e', 'lam'])
al, be, e, lam = K.gens
one = K(1)
ctx = K.ctx


def poly_pos(p):
    """all coefficients > 0 (p a fmpq_mpoly)"""
    return all(c > 0 for c in p.coeffs())


def bernstein_pos(p, ilam=3):
    """write p (fmpq_mpoly) as sum_k p_k lam^k; Bernstein coefficients of degree d = deg_lam p in lam;
    return (all coefficients of every Bernstein coefficient >= 0, number of zero Bernstein coefficients, min #terms)."""
    d = p.degrees()[ilam]
    parts = [ctx.from_dict({}) for _ in range(d + 1)]
    for exps, cf in zip(p.monoms(), p.coeffs()):
        k = exps[ilam]
        ex = list(exps); ex[ilam] = 0
        parts[k] = parts[k] + ctx.from_dict({tuple(ex): cf})
    ok = True; nz = 0
    from math import comb
    for j in range(d + 1):
        b = ctx.from_dict({})
        for kk in range(j + 1):
            b = b + parts[kk] * flint.fmpq(comb(j, kk), comb(d, kk))
        if b.is_zero():
            nz += 1
            if j == 0 or j == d:
                ok = False   # positivity on the closed interval needs the end coefficients nonzero
            continue
        if not all(c >= 0 for c in b.coeffs()):
            ok = False
    return ok, nz, d


def factors_pos(q):
    c, fac = q.factor()
    okf = c > 0
    out = []
    for f, ex in fac:
        if poly_pos(f):
            out.append((f.total_degree(), len(f), ex, 'pos'))
        else:
            # allow factors in lam that are positive on [0,1] via Bernstein
            okb, nz, d = bernstein_pos(f)
            if okb and nz <= d:
                out.append((f.total_degree(), len(f), ex, 'bernstein'))
            elif ex % 2 == 0:
                out.append((f.total_degree(), len(f), ex, 'even power'))
            else:
                okf = False
                out.append((f.total_degree(), len(f), ex, 'NOT POSITIVE'))
    return okf, out


ga = al / (2 * (be + e)); gb = be / (2 * (al + e))
kap = (1 + ga) * (1 + gb)
G0 = 1 + (ga + gb) / 2
ap = (al + be + 2 * e) / (be + 2 * e)
bp = (al + be + 2 * e) / (al + 2 * e)
K1 = (1 + ap) / 2; K2 = (1 + bp) / 2
Z0 = one / (al + be + 2 * e); Z1 = one / (al + be)
Z = Z0 + lam * (Z1 - Z0)
ds = [al + be, al, be]


def h(d, z):
    return z / (1 + d * z)


Frng = []
for d in ds:
    hi = e * h(d, Z)
    lo = e * ((1 - lam) * h(d, Z0) + lam * h(d, Z1))
    Frng.append((lo, hi))
G3 = 1 + e * Z / 2
ok = True
# c < 2 at every F0 in range: 2 - c >= 2 - kap/(2 G0)
cmax = kap / (2 * G0)
r = 2 - cmax
okc = poly_pos(r.n) and poly_pos(r.d)
print(f'c < 2: 2 - kap/(2 G0) has numerator and denominator with positive coefficients: {okc}')
ok &= okc
for i0, i1, i2 in itertools.product((0, 1), repeat=3):
    c = kap / (2 * (G0 + Frng[0][i0] / 2))
    A = ap / (K1 + Frng[1][i1] / 2)
    B = bp / (K2 + Frng[2][i2] / 2)
    g = 4 - c * c - G3 * (A * A + B * B - c * A * B)
    okN, nz, d = bernstein_pos(g.n)
    okD, facs = factors_pos(g.d)
    tag = ''.join('lh'[i] for i in (i0, i1, i2))
    print(f'F0,F1,F2 = {tag}: numerator {len(g.n)} terms, total degree {g.n.total_degree()}, lam-degree {d}; '
          f'Bernstein coefficients nonnegative: {okN} (zero ones: {nz}); denominator factors {facs}; ok {okD} '
          f'[{time.time() - t0:.1f}s]', flush=True)
    ok &= okN and okD and nz <= d
print('PASS' if ok else 'FAIL')
