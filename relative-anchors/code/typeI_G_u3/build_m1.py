"""u = 3, m = 1: polynomial data for G > 0 (closed form from closed_m1.py, z-coordinates).

G_z = Nt/den with Nt polynomial and one common denominator den (checked); integer coefficients after multiplying Nt and
den by one positive integer L.  G_z > 0 on the open orthant is equivalent to G > 0 (G_z = Dg G_zeta Dg,
Dg = diag(1, w0 w1, w0 w2, w0 w3) invertible).  The leading principal minors of G_z are minor_j = m_j / den^j with
m_j = det Nt[:j, :j].  Write den = cst * prod_i f_i^{E_i} (irreducible factors).  We compute every m_j as
prod_i f_i^{x_ij} * c_j with the largest exponents x_ij (exact divisions, checked by divmod), by Laplace expansion:
  - the 2 x 2 minors M(01|ab) of rows {0, 1} and M(23|ab) of rows {2, 3} are stripped of the factors f_i;
  - m_3 = sum_j +- Nt_2j M(01|..) and m_4 = sum +- M(01|ab) M(23|cd) (generalised Laplace expansion), with the
    common factors pulled out before summing.
The expansions are checked against direct determinants of Nt at a rational point.  Output: build_m1.pkl."""
import pickle, sys, time, flint
from math import lcm
from fractions import Fraction as Fr

u = 3
t0 = time.time()
D = pickle.load(open(f'closed_m1_u{u}.pkl', 'rb'))
names = tuple(D['names'])
ctxZ = flint.fmpz_mpoly_ctx.get(names, 'lex')
G = D['G']
d = len(G)
L = 1
for row in G:
    for (n_, d_) in row:
        for c in list(n_.values()) + list(d_.values()):
            L = lcm(L, int(c.q))


def toZ(dct):
    out = {}
    for m, c in dct.items():
        v = c * L
        assert v.q == 1
        out[m] = int(v.p)
    return ctxZ.from_dict(out)


Nt = [[toZ(G[i][j][0]) for j in range(d)] for i in range(d)]
den = toZ(G[0][0][1])
assert all(toZ(G[i][j][1]) == den for i in range(d) for j in range(d))
assert all(Nt[i][j] == Nt[j][i] for i in range(d) for j in range(d))
print(f'loaded, L = {L} ({time.time() - t0:.1f} s)', flush=True)


def divout(p, q):
    e = 0
    while True:
        qq, r = divmod(p, q)
        if r != 0:
            return p, e
        p = qq
        e += 1


cst, fac = den.factor()
FACS = [f for f, e in fac]
print('den = ' + str(cst) + ' * ' + ' * '.join(f'[deg {f.total_degree()}, {len(f)} terms]^{e}' for f, e in fac), flush=True)


def strip_facs(p):
    """p = prod f_i^{e_i} * core; returns (core, [e_i])."""
    ex = []
    for f in FACS:
        p, e = divout(p, f)
        ex.append(e)
    return p, ex


def common(items):
    """items: list of (core, ex).  Returns (g_ex, [core * prod f^(ex - g_ex)])."""
    g = [min(it[1][i] for it in items) for i in range(len(FACS))]
    out = []
    for c, ex in items:
        p = c
        for f, e, gi in zip(FACS, ex, g):
            if e > gi:
                p = p * f ** (e - gi)
        out.append(p)
    return g, out


def minor2(r1, r2, a, b):
    return Nt[r1][a] * Nt[r2][b] - Nt[r1][b] * Nt[r2][a]


pairs = [(a, b) for a in range(4) for b in range(a + 1, 4)]
S01 = {ab: strip_facs(minor2(0, 1, *ab)) for ab in pairs}
S23 = {ab: strip_facs(minor2(2, 3, *ab)) for ab in pairs}
print(f'2 x 2 minors stripped ({time.time() - t0:.1f} s): ' + ', '.join(f'{ab}: {S01[ab][1]}' for ab in pairs), flush=True)
# m_1, m_2 as (core, exponents)
m1 = strip_facs(Nt[0][0])
m2 = S01[(0, 1)]
# m_3 = Nt20 M01(1,2) - Nt21 M01(0,2) + Nt22 M01(0,1)
g3, (A12, A02, A01) = common([S01[(1, 2)], S01[(0, 2)], S01[(0, 1)]])
c3, e3 = strip_facs(Nt[2][0] * A12 - Nt[2][1] * A02 + Nt[2][2] * A01)
m3 = (c3, [x + y for x, y in zip(g3, e3)])
# m_4 = sum over a < b of (-1)^(a + b + 1) M01(a,b) M23(complement)
gA, LA = common([S01[ab] for ab in pairs])
gB, LB = common([S23[tuple(x for x in range(4) if x not in ab)] for ab in pairs])
tot = 0 * Nt[0][0]
for k, ab in enumerate(pairs):
    tot = tot + (-1) ** (ab[0] + ab[1] + 1) * LA[k] * LB[k]
c4, e4 = strip_facs(tot)
m4 = (c4, [x + y + z for x, y, z in zip(gA, gB, e4)])
print(f'expansions done ({time.time() - t0:.1f} s)', flush=True)


def full(mm):
    c, ex = mm
    p = c
    for f, e in zip(FACS, ex):
        p = p * f ** e
    return p


import random
random.seed(3)
pt = [random.randint(1, 30) for _ in range(len(names))]
Nv = [[int(Nt[i][j](*pt)) for j in range(d)] for i in range(d)]


def detF(A):
    A = [[Fr(x) for x in row] for row in A]
    n = len(A)
    det_ = Fr(1)
    for p in range(n):
        piv = next(i for i in range(p, n) if A[i][p] != 0)
        if piv != p:
            A[p], A[piv] = A[piv], A[p]
            det_ = -det_
        det_ *= A[p][p]
        for i in range(p + 1, n):
            f = A[i][p] / A[p][p]
            for j in range(p, n):
                A[i][j] -= f * A[p][j]
    return det_


okc = True
for j, mm in enumerate((m1, m2, m3, m4)):
    lhs = detF([row[:j + 1] for row in Nv[:j + 1]])
    rhs = Fr(int(full(mm)(*pt)))
    okc &= lhs == rhs
    print(f'  check m_{j + 1} = det Nt[:{j + 1}, :{j + 1}] at a rational point: {lhs == rhs}', flush=True)
cores = []
for j, mm in enumerate((m1, m2, m3, m4)):
    cores.append(mm)
    print(f'minor {j + 1}: m_{j + 1} = prod(den factors^{mm[1]}) * core, core degree {mm[0].total_degree()}, {len(mm[0])} terms', flush=True)
for j, (cj, e) in enumerate(cores):
    cc, fc = cj.factor()
    print(f'core {j + 1} = {cc} * ' + ' * '.join(f'[deg {f.total_degree()}, {len(f)} terms]^{e}' for f, e in fc), flush=True)
res = dict(L=L, names=names, den={m: int(v) for m, v in den.to_dict().items()},
           cores=[({m: int(v) for m, v in cj.to_dict().items()}, e) for cj, e in cores],
           den_factors=[({m: int(v) for m, v in f.to_dict().items()}, e) for f, e in fac], den_const=int(cst))
pickle.dump(res, open('build_m1.pkl', 'wb'))
print(f'time {time.time() - t0:.1f} s')
print('PASS' if okc else 'FAIL')
