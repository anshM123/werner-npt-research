"""Step 2 of the corner proof: from the 15 node determinants det Q7(q0, t, nu; tau_l), tau_l = 2..16, recover
A * det Q7 as a polynomial in (q0, t, nu, tau) by exact Lagrange interpolation in tau (degree <= 14 in tau, by the
structure lemma checked in check_structure.py), split off the explicit factor
    L = tau^3 * i^2 * (i + q0 + t + 1)^2,   i = (nu - tau)/2,
and the tau-free content g(q0, t, nu), and save the core polynomial P (integer coefficients) in core.pkl.
det Q7 = c * L * g * P / A  with c a rational constant."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
import pickle, time, sys
from fractions import Fraction
import flint

T0 = time.time()
V3 = ('q0', 't', 'nu')
V4 = ('q0', 't', 'nu', 'tau')
ctx3 = flint.fmpq_mpoly_ctx.get(V3, 'lex')
ctx4 = flint.fmpq_mpoly_ctx.get(V4, 'lex')
taus = list(range(2, 17))


def load(tau):
    with open(f'nodes/det_tau{tau}.pkl', 'rb') as f:
        D = pickle.load(f)
    assert tuple(D['vars']) == V3 and D['tau'] == tau
    n = ctx3.from_dict({e: flint.fmpq(p, q) for e, (p, q) in D['n'].items()})
    d = ctx3.from_dict({e: flint.fmpq(p, q) for e, (p, q) in D['d'].items()})
    return n, d


nodes = {tau: load(tau) for tau in taus}
A = nodes[taus[0]][1]
same = all(nodes[tau][1] == A for tau in taus)
print('node denominators all equal:', same, '; A terms', len(A), 'total degree', A.total_degree(), flush=True)
if not same:
    for tau in taus[1:]:
        d = nodes[tau][1]
        g = A.gcd(d)
        A = A * (d / g)
    print('  lcm A terms', len(A), 'total degree', A.total_degree(), flush=True)
N_l = {}
for tau in taus:
    n, d = nodes[tau]
    q = A / d                      # exact division (A is a common multiple)
    assert q * d == A
    N_l[tau] = n * q
print('scaled node numerators', f'{time.time()-T0:.1f}s', flush=True)


# Lagrange basis coefficients (exact rationals)
def polymul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for x, u in enumerate(a):
        for y, v in enumerate(b):
            out[x + y] += u * v
    return out


nb = len(taus)
coef = {d: ctx3.from_dict({}) for d in range(nb)}
for l, tl in enumerate(taus):
    basis = [Fraction(1)]
    for mm, tm in enumerate(taus):
        if mm != l:
            basis = polymul(basis, [Fraction(-tm, tl - tm), Fraction(1, tl - tm)])
    for d in range(nb):
        if basis[d] != 0:
            coef[d] = coef[d] + N_l[tl] * flint.fmpq(basis[d].numerator, basis[d].denominator)
degs = [d for d in range(nb) if not coef[d].is_zero()]
print('tau-degree of A*det:', max(degs), '; nonzero tau-powers', degs, f'{time.time()-T0:.1f}s', flush=True)
# 4-variable polynomial N = A * det
Nd = {}
for d in degs:
    for e, c in coef[d].to_dict().items():
        Nd[e + (d,)] = c
N = ctx4.from_dict(Nd)
print('N terms', len(N), 'degrees', N.degrees(), 'total', N.total_degree(), flush=True)
q0, t, nu, tau = ctx4.gens()
i_ = (nu - tau) / 2
L = tau ** 3 * i_ ** 2 * (i_ + q0 + t + 1) ** 2
R = N / L
assert R * L == N
print('N divisible by L = tau^3 i^2 (i+m+1)^2: True; R terms', len(R), 'degrees', R.degrees(), flush=True)
# tau-free content: gcd of the tau-coefficients of R
Rc = {}
for e, c in R.to_dict().items():
    Rc.setdefault(e[3], {})[e[:3]] = c
Rcoef = {d: ctx3.from_dict(v) for d, v in Rc.items()}
g = None
for d in sorted(Rcoef, key=lambda d: len(Rcoef[d])):
    g = Rcoef[d] if g is None else g.gcd(Rcoef[d])
    if g.total_degree() == 0:
        break
print('tau-content g =', g, f'({time.time()-T0:.1f}s)', flush=True)
gd = {e + (0,): c for e, c in g.to_dict().items()}
g4 = ctx4.from_dict(gd)
P = R / g4
assert P * g4 == R
# integer normalisation: P = (cont/den) * Pz with Pz a primitive integer polynomial
from math import gcd
Pd = P.to_dict()
den = 1
for c in Pd.values():
    den = den * int(c.q) // gcd(den, int(c.q))
Pz = {e: int(c.p) * (den // int(c.q)) for e, c in Pd.items()}
cont = 0
for v in Pz.values():
    cont = gcd(cont, abs(v))
Pz = {e: v // cont for e, v in Pz.items()}
print('core P: terms', len(Pz), 'degrees', P.degrees(), 'total', P.total_degree(), f'{time.time()-T0:.1f}s', flush=True)
with open('core.pkl', 'wb') as f:
    pickle.dump({'vars': V4, 'P': Pz, 'g': {e: (int(c.p), int(c.q)) for e, c in g.to_dict().items()},
                 'A': {e: (int(c.p), int(c.q)) for e, c in A.to_dict().items()},
                 'N': {e: (int(c.p), int(c.q)) for e, c in N.to_dict().items()}}, f)
print('saved core.pkl', f'{time.time()-T0:.1f}s')
