"""Certificate for the blocks of (R_inf)_2 with |T| = 3 (no letter 1), every k >= 4 (NOTES.md, Section 5, Theorem D2).

Homogeneous coordinates: a1, a2, a3 = x_t (t in T), e = sum of the other x_u (a1 + a2 + a3 + e = |s| = 1 after
scaling; every quantity below is homogeneous of degree 0).  xT = a1 + a2 + a3, Z = F_T/e in [Z0, Z1],
Z0 = 1/(xT + 2e), Z1 = 1/xT, Z = Z0 + lam (Z1 - Z0).  For the pair {i, j} with third index k_ and for i:
    c_ij = a'_i a'_j / (2 (G0_ij + F_ij/2)),  a'_i = 1 + g(a_i),  g(u) = u/(2 - 2u - a_k_),  G0_ij = 1 + (g(a_i) + g(a_j))/2,
    F_ij in [e l_d(Z), e h_d(Z)] with d = xT - a_k_;
    A_i = (1 + g'_i)/(1 + (g'_i + F'_i)/2),  g'_i = a_i/(2 - a_i - xT),  F'_i in [e l_d(Z), e h_d(Z)] with d = a_i;
    G_T = 1 + e Z/2,  h_d(z) = z/(1 + d z),  l_d the chord of h_d on [Z0, Z1].
Checked (exact, python-flint): (1) at each of the 8 vertices of the c-box, the leading minors 4 - c12^2 and det C of
C = [[2, c12, c13], [c12, 2, c23], [c13, c23, 2]] are positive; (2) at each of the 64 vertices of the (c, A)-box,
P = 2 det C - G_T a^T adj(C) a is positive.  Positivity of a rational function is certified as in cert_t2.py:
numerator in the Bernstein basis of lam whose coefficients are polynomials in (a1, a2, a3, e) with only nonnegative
coefficients, the first and the last of them nonzero (so the numerator is positive for a1, a2, a3, e > 0 and every
lam in [0, 1]); denominator a product of factors with only positive coefficients.
Usage: python cert_t3j2.py"""
import itertools, time
from math import comb
import flint
from flintfield import Field

t0 = time.time()
K = Field(['a1', 'a2', 'a3', 'e', 'lam'])
a1, a2, a3, e, lam = K.gens
ctx = K.ctx
ILAM = 4
A = [a1, a2, a3]
Sg = a1 + a2 + a3 + e
xT = a1 + a2 + a3
Z0 = K(1) / (xT + 2 * e)
Z1 = K(1) / xT
Z = Z0 + lam * (Z1 - Z0)
GT = 1 + e * Z / 2


def h(d, z):
    return z / (1 + d * z)


def Frange(d):
    return (e * ((1 - lam) * h(d, Z0) + lam * h(d, Z1)), e * h(d, Z))


def bernstein_pos(p):
    d = p.degrees()[ILAM]
    parts = [ctx.from_dict({}) for _ in range(d + 1)]
    for exps, cf in p.to_dict().items():
        kk = exps[ILAM]
        ex = list(exps); ex[ILAM] = 0
        parts[kk] = parts[kk] + ctx.from_dict({tuple(ex): cf})
    ok = True; nz = 0
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


def den_pos(q):
    c, fac = q.factor()
    return c > 0 and all(all(x > 0 for x in f.coeffs()) for f, _ in fac), [(f.total_degree(), len(f), ex) for f, ex in fac]


def certify(r, tag):
    okn, nz, d = bernstein_pos(r.n)
    okd, facs = den_pos(r.d)
    return okn and okd, f'{tag}: numerator {len(r.n)} terms, degree {r.n.total_degree()}, lam-degree {d}, ' \
                        f'Bernstein ok {okn} (zero {nz}); denominator factors {facs} ok {okd}'


pairs = [(0, 1, 2), (0, 2, 1), (1, 2, 0)]   # (i, j, third)
crng = []
for i, j, kk in pairs:
    g = lambda u: u / (2 * Sg - 2 * u - A[kk])
    api, apj = 1 + g(A[i]), 1 + g(A[j])
    G0 = 1 + (g(A[i]) + g(A[j])) / 2
    crng.append([api * apj / (2 * (G0 + F / 2)) for F in Frange(xT - A[kk])])
arng = []
for i in range(3):
    gp = A[i] / (2 * Sg - A[i] - xT)
    arng.append([(1 + gp) / (1 + (gp + F) / 2) for F in Frange(A[i])])
ok = True
# (1) C > 0 at the 8 c-vertices
for cv in itertools.product((0, 1), repeat=3):
    c12, c13, c23 = crng[0][cv[0]], crng[1][cv[1]], crng[2][cv[2]]
    m2 = 4 - c12 * c12
    det = 8 - 2 * (c12 * c12 + c13 * c13 + c23 * c23) + 2 * c12 * c13 * c23
    o1, s1 = certify(m2, f'C minor 2 at {cv}')
    o2, s2 = certify(det, f'det C at {cv}')
    print(s1, flush=True); print(s2, f'[{time.time() - t0:.1f}s]', flush=True)
    ok &= o1 and o2
# (2) P > 0 at the 64 vertices
cnt = 0
for cv in itertools.product((0, 1), repeat=3):
    c12, c13, c23 = crng[0][cv[0]], crng[1][cv[1]], crng[2][cv[2]]
    det = 8 - 2 * (c12 * c12 + c13 * c13 + c23 * c23) + 2 * c12 * c13 * c23
    adj = [[4 - c23 * c23, c13 * c23 - 2 * c12, c12 * c23 - 2 * c13],
           [c13 * c23 - 2 * c12, 4 - c13 * c13, c12 * c13 - 2 * c23],
           [c12 * c23 - 2 * c13, c12 * c13 - 2 * c23, 4 - c12 * c12]]
    for av in itertools.product((0, 1), repeat=3):
        a = [arng[i][av[i]] for i in range(3)]
        q = sum((a[i] * adj[i][j] * a[j] for i in range(3) for j in range(3)), K(0))
        P = 2 * det - GT * q
        o, s = certify(P, f'P at c {cv}, A {av}')
        cnt += 1
        print(s, f'[{time.time() - t0:.1f}s]', flush=True)
        ok &= o
print(f'{cnt} vertex inequalities checked')
print('PASS' if ok else 'FAIL')
