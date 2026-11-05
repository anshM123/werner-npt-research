"""Checks behind the interpolation step of the corner proof (all exact, python-flint):
 (1) window exactness of the generic calculus on the corner, and Q8, G8 of SymF2 over QQ(q0, t, i, tau) agree with
     the numeric calculus Num2 at 3^4 = 81 sample points of the corner;
 (2) STRUCTURE LEMMA: written in (q0, t, nu, tau), nu = 2i + tau, every entry of Q7 has a tau-free denominator and a
     numerator of degree <= 2 in tau;  hence A * det Q7 is a polynomial of degree <= 14 in tau for a tau-free A;
 (3) the node matrices of compute_nodes.py (SymNu at tau = 2..16) are the specialisations of the same Q7;
 (4) every irreducible denominator factor ('atom') of Q7 is positive on the corner (one-signed coefficients after
     q0 = 2 + x, t = 6 + y, nu = 6 + v), and the common denominator A of the node determinants is a product of atoms."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
import itertools, time, pickle, sys
from math import comb
from fractions import Fraction
import flint
from orbit2 import SymF2, Num2, reduced2
from nucalc import SymNu, form_in_basis, Q7NAMES, V8c

T0 = time.time()
S0, S1, S2, S3 = (int(v) for v in sys.argv[1:5]) if len(sys.argv) > 4 else (2, 6, 2, 2)
print(f'region: q0 >= {S0}, t >= {S1}, i >= {S2}, tau >= {S3}', flush=True)
assert S0 >= 2 and S1 >= 4 and S2 >= 2 and S3 >= 1      # inside the generic family
ok_all = True
# ---------------------------------------------------------------- (1) window + sanity on the corner
be = SymF2(None, None, None, None)
V8n = list(V8c)
Q8, G8 = form_in_basis(be, V8n)
dmin = min(d for d, dp in be.seen); dpmin = min(dp for d, dp in be.seen)
win = dmin + S0 >= 2 and dpmin + S1 >= 2        # every node genuine (offsets -2..2 valid) on the region
print(f'(1) level offsets: delta >= {dmin}, delta\' >= {dpmin}; window exact on the region: {win}', flush=True)
ok_all &= win


def evalF(x, vals):
    if x == 0:
        return Fraction(0)
    r = x.n(*vals) / x.d(*vals)
    return Fraction(int(r.p), int(r.q))


nsamp = 0
for (q0, t, i, tau) in itertools.product((S0, S0 + 1, S0 + 3), (S1, S1 + 1, S1 + 3), (S2, S2 + 1, S2 + 3), (S3, S3 + 1, S3 + 3)):
    m = q0 + t; k = m + 2 * i + tau; j = q0 + i + tau
    names, G, Q = reduced2(Num2(k, j, tau, i), 'sym')
    N0 = 2 ** i * comb(m, q0)
    idx = {nm: names.index(nm) for nm in names}
    for a, na in enumerate(V8n):
        for b, nb in enumerate(V8n):
            qn = sum(ca * cb * Q[idx[x]][idx[y]] for x, ca in V8c[na].items() for y, cb in V8c[nb].items())
            gn = sum(ca * cb * G[idx[x]][idx[y]] for x, ca in V8c[na].items() for y, cb in V8c[nb].items())
            if evalF(Q8[a][b], (q0, t, i, tau)) * N0 != qn or evalF(G8[a][b], (q0, t, i, tau)) * N0 != gn:
                ok_all = False
                print('   SANITY FAILURE at', (q0, t, i, tau), na, nb)
    nsamp += 1
print(f'(1) SymF2 = Num2 (Q8 and Gram, all 64 entries) at {nsamp} sample points of the region: {ok_all}  ({time.time()-T0:.1f}s)', flush=True)
# ---------------------------------------------------------------- (2) structure lemma
Q7, G7 = form_in_basis(be, Q7NAMES)
ctx4 = flint.fmpq_mpoly_ctx.get(('q0', 't', 'nu', 'tau'), 'lex')
Q0, Tt, NU, TAU = ctx4.gens()
I_ = (NU - TAU) / 2
sub = (Q0, Tt, I_, TAU)


def to_nu(x):
    """entry over QQ(q0, t, i, tau) -> reduced (num, den) over QQ[q0, t, nu, tau]"""
    if x == 0:
        return ctx4.from_dict({}), ctx4.from_dict({(0, 0, 0, 0): 1})
    n = x.n.compose(*sub, ctx=ctx4); d = x.d.compose(*sub, ctx=ctx4)
    g = n.gcd(d)
    return n / g, d / g


E7 = [[to_nu(Q7[a][b]) for b in range(7)] for a in range(7)]
struct = all(E7[a][b][1].degrees()[3] == 0 and E7[a][b][0].degrees()[3] <= 2 for a in range(7) for b in range(7))
print(f'(2) structure lemma (tau-free denominators, tau-degree <= 2) for all 49 entries of Q7: {struct}', flush=True)
ok_all &= struct
# ---------------------------------------------------------------- (3) nodes are specialisations
ctx3 = flint.fmpq_mpoly_ctx.get(('q0', 't', 'nu'), 'lex')
g3 = ctx3.gens()
ok3 = True
for tau in range(2, 17):
    bn = SymNu(tau)
    Qn, Gn = form_in_basis(bn, Q7NAMES)
    for a in range(7):
        for b in range(7):
            n, d = E7[a][b]
            ns = n.compose(g3[0], g3[1], g3[2], ctx3.from_dict({(0, 0, 0): tau}), ctx=ctx3)
            ds = d.compose(g3[0], g3[1], g3[2], ctx3.from_dict({(0, 0, 0): tau}), ctx=ctx3)
            x = Qn[a][b]
            if x.n * ds != ns * x.d:
                ok3 = False
print(f'(3) SymNu(tau) = Q7 at tau for tau = 2..16 (all entries): {ok3}  ({time.time()-T0:.1f}s)', flush=True)
ok_all &= ok3
# ---------------------------------------------------------------- (4) atoms
atoms = []
for a in range(7):
    for b in range(7):
        d = E7[a][b][1]
        if d.total_degree() == 0:
            continue
        for f, e in d.factor()[1]:
            if all(f != h for h in atoms):
                atoms.append(f)
print(f'(4) {len(atoms)} distinct atoms:', [str(f) for f in atoms], flush=True)
x_, y_, v_ = g3
ok4 = True
atoms3 = []
for f in atoms:
    # f is tau-free: move to QQ[q0, t, nu]
    f3 = ctx3.from_dict({e[:3]: c for e, c in f.to_dict().items()})
    atoms3.append(f3)
    s = f3.compose(x_ + S0, y_ + S1, v_ + (2 * S2 + S3))
    co = list(s.coeffs())
    c0 = s.to_dict().get((0, 0, 0), 0)
    one = (all(c > 0 for c in co) and c0 > 0) or (all(c < 0 for c in co) and c0 < 0)
    ok4 &= one
print(f'(4) every atom one-signed with nonzero constant after the shift (no zero on the real region): {ok4}', flush=True)
ok_all &= ok4
try:
    with open('nodes/det_tau2.pkl', 'rb') as f:
        D = pickle.load(f)
    A = ctx3.from_dict({e: flint.fmpq(p, q) for e, (p, q) in D['d'].items()})
    rem = A
    expo = []
    for f3 in atoms3:
        e = 0
        while True:
            qq, rr = divmod(rem, f3)
            if rr.is_zero():
                rem = qq; e += 1
            else:
                break
        expo.append(e)
    okA = rem.total_degree() == 0 and not rem.is_zero()
    print(f'(4) common denominator A = const * prod atoms^e, e = {expo}, const = {rem}: {okA}', flush=True)
    ok_all &= okA
except FileNotFoundError:
    print('(4) nodes/det_tau2.pkl not found: A not checked')
    ok_all = False
print('STRUCTURE CHECKS PASS' if ok_all else 'STRUCTURE CHECKS FAIL', f'{time.time()-T0:.0f}s')
