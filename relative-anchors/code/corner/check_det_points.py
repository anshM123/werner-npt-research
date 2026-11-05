"""Cross-checks of the interpolated determinant (core.pkl) against direct exact determinants, and positive
definiteness of Q7 at one point of the corner:
 (a) at 40 random rational points (q0, t, i, tau) and at 30 integer points of the corner, det Q7 (entries of the
     4-variable SymF2 matrix evaluated exactly, then an exact 7 x 7 determinant) equals N/A, and equals
     c * L * g * P / A with ONE rational constant c (L = tau^3 i^2 (i + q0 + t + 1)^2);
 (b) at the corner point (q0, t, i, tau) = (2, 6, 2, 2) (k = 14), Q7 is positive definite: exact LDL of the matrix
     built by the numeric calculus Num2 (independent of the symbolic backend), all 7 pivots > 0;
 (c) at six larger integer points (k up to 44) the determinant of the Num2 matrix, scaled back by 2^i C(m, q0),
     equals N/A and is positive.
Usage: python check_det_points.py [q0 t i tau] (the base point; default 2 6 2 2, the corner's vertex)."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
import pickle, random, time, sys
from fractions import Fraction
from math import comb
import flint
from orbit2 import SymF2, Num2, reduced2
from nucalc import form_in_basis, Q7NAMES, V8c

T0 = time.time()
PT = tuple(int(v) for v in sys.argv[1:5]) if len(sys.argv) > 4 else (2, 6, 2, 2)
with open('core.pkl', 'rb') as f:
    D = pickle.load(f)
ctx4 = flint.fmpq_mpoly_ctx.get(('q0', 't', 'nu', 'tau'), 'lex')
ctx3 = flint.fmpq_mpoly_ctx.get(('q0', 't', 'nu'), 'lex')
N = ctx4.from_dict({e: flint.fmpq(p, q) for e, (p, q) in D['N'].items()})
A = ctx3.from_dict({e: flint.fmpq(p, q) for e, (p, q) in D['A'].items()})
g = ctx3.from_dict({e: flint.fmpq(p, q) for e, (p, q) in D['g'].items()})
P = ctx4.from_dict({e: flint.fmpq(v) for e, v in D['P'].items()})
be = SymF2(None, None, None, None)
Q7, G7 = form_in_basis(be, Q7NAMES)


def ev(x, pt):
    if x == 0:
        return flint.fmpq(0)
    return x.n(*pt) / x.d(*pt)


def check_point(pt):
    q0, t, i, tau = pt
    nu = 2 * i + tau
    M = flint.fmpq_mat(7, 7, [ev(Q7[a][b], pt) for a in range(7) for b in range(7)])
    det = M.det()
    rhs = N(q0, t, nu, tau) / A(q0, t, nu)
    L = tau ** 3 * i ** 2 * (i + q0 + t + 1) ** 2
    fac = L * g(q0, t, nu) * P(q0, t, nu, tau) / A(q0, t, nu) if A(q0, t, nu) != 0 else None
    return det, rhs, fac


random.seed(1)
ok = True
c = None
pts = [tuple(flint.fmpq(random.randint(-50, 50), random.randint(1, 9)) for _ in range(4)) for _ in range(40)]
pts += [(random.randint(2, 30), random.randint(6, 30), random.randint(2, 20), random.randint(2, 40)) for _ in range(30)]
nchk = 0
for pt in pts:
    pt = tuple(flint.fmpq(x) for x in pt)
    try:
        det, rhs, fac = check_point(pt)
    except ZeroDivisionError:
        continue                      # a random point on a pole: skip
    nchk += 1
    if det != rhs:
        ok = False; print('  MISMATCH det vs N/A at', pt)
    if fac is not None and fac != 0:
        r = det / fac
        if c is None:
            c = r
        elif r != c:
            ok = False; print('  MISMATCH of the constant at', pt)
print(f'(a) det Q7 = N/A = c L g P / A at {nchk} points (c = {c}): {ok}  ({time.time()-T0:.1f}s)', flush=True)
# (b) PD at the point PT via Num2
q0, t, i, tau = PT
m = q0 + t; k = m + 2 * i + tau; j = q0 + i + tau
names, G, Q = reduced2(Num2(k, j, tau, i), 'sym')
idx = {nm: names.index(nm) for nm in names}
M = [[sum(ca * cb * Q[idx[x]][idx[y]] for x, ca in V8c[na].items() for y, cb in V8c[nb].items()) for nb in Q7NAMES]
     for na in Q7NAMES]
piv = []
M = [row[:] for row in M]
for x in range(7):
    p = M[x][x]
    piv.append(p)
    if p <= 0:
        break
    for r in range(x + 1, 7):
        f = M[r][x] / p
        for s in range(x + 1, 7):
            M[r][s] -= f * M[x][s]
pd = len(piv) == 7 and all(p > 0 for p in piv)
print(f'(b) Q7 at (q0, t, i, tau) = {PT}, k = {k}: pivots {[float(p) for p in piv]} -> positive definite: {pd}')
ok &= pd
# also the symbolic matrix at that point has the same determinant sign
det, rhs, fac = check_point(tuple(flint.fmpq(x) for x in PT))
print(f'    det there (symbolic matrix) = {float(det):.4e}, N/A = {float(rhs):.4e}')
# (c) independent of the symbolic backend: det of the Num2 matrix (scaled back by 2^i C(m, q0)) = N/A at larger
#     integer points of the region
okc = True
extra = [(PT[0] + a, PT[1] + b, PT[2] + c_, PT[3] + d) for (a, b, c_, d) in
         [(5, 5, 4, 7), (1, 14, 2, 11), (10, 0, 7, 1), (0, 0, 0, 20), (3, 2, 12, 0), (8, 9, 1, 2)]]
for (q0, t, i, tau) in extra:
    m = q0 + t; k = m + 2 * i + tau; j = q0 + i + tau
    names, G, Q = reduced2(Num2(k, j, tau, i), 'sym')
    idx = {nm: names.index(nm) for nm in names}
    Mn = [[sum(ca * cb * Q[idx[x]][idx[y]] for x, ca in V8c[na].items() for y, cb in V8c[nb].items()) for nb in Q7NAMES]
          for na in Q7NAMES]
    N0 = 2 ** i * comb(m, q0)
    dn = flint.fmpq_mat(7, 7, [flint.fmpq(x.numerator, x.denominator) for row in Mn for x in row]).det() / flint.fmpq(N0) ** 7
    nu = 2 * i + tau
    rhs = N(q0, t, nu, tau) / A(q0, t, nu)
    okc &= (dn == rhs) and dn > 0
    print(f'    (c) k = {k}, (q0, t, i, tau) = {(q0, t, i, tau)}: det(Num2) = N/A: {dn == rhs}, positive: {dn > 0}')
print(f'(c) Num2 determinants agree with N/A at {len(extra)} larger points: {okc}')
ok &= okc
print('POINT CHECKS PASS' if ok else 'POINT CHECKS FAIL', f'{time.time()-T0:.0f}s')
