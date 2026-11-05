"""Step 3 of the corner proof: the core polynomial P(q0, t, nu, tau) (core.pkl) has all coefficients of one sign
after the shift to the corner's box variables
    q0 = 2 + x,  t = 6 + y,  i = 2 + z,  tau = 2 + w   (so nu = 2i + tau = 6 + 2z + w),   x, y, z, w >= 0,
and a nonzero constant term; hence P has no zero on the real corner [2, oo) x [6, oo) x [2, oo) x [2, oo).
Done with bivariate operations only (exact, python-flint fmpz_mpoly):
  (1) P = sum_{c,d} nu^c tau^d p_cd(q0, t); shift each p_cd(S0 + x, S1 + y);
  (2) regroup as sum_{a,b} x^a y^b r_ab(nu, tau) and expand each r_ab(2 S2 + S3 + 2z + w, S3 + w).
Usage: python shift_check.py [core.pkl] [q0min tmin imin taumin]   (default 2 6 2 2 = the corner)"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
import pickle, sys, time
import flint

T0 = time.time()
fn = sys.argv[1] if len(sys.argv) > 1 else 'core.pkl'
S0, S1, S2, S3 = (int(v) for v in sys.argv[2:6]) if len(sys.argv) > 5 else (2, 6, 2, 2)
SNU = 2 * S2 + S3
print(f'shift: q0 = {S0} + x, t = {S1} + y, i = {S2} + z, tau = {S3} + w  (nu = {SNU} + 2z + w)', flush=True)
with open(fn, 'rb') as f:
    D = pickle.load(f)
assert tuple(D['vars']) == ('q0', 't', 'nu', 'tau')
Pd = D['P']
print('P terms', len(Pd), 'max exponents', [max(e[k] for e in Pd) for k in range(4)],
      'total degree', max(sum(e) for e in Pd), flush=True)
ctxqt = flint.fmpz_mpoly_ctx.get(('q0', 't'), 'lex')
X, Y = ctxqt.gens()
s1 = {}
for e, c in Pd.items():
    s1.setdefault((e[2], e[3]), {})[(e[0], e[1])] = c
s2 = {}
for (c_, d_), dct in s1.items():
    p = ctxqt.from_dict(dct).compose(X + S0, Y + S1)
    for (a, b), cc in p.to_dict().items():
        s2.setdefault((a, b), {})[(c_, d_)] = cc
print(f'(1) shifted (q0, t) in {len(s1)} slices; {len(s2)} (x, y)-slices ({time.time()-T0:.1f}s)', flush=True)
ctxnt = flint.fmpz_mpoly_ctx.get(('nu', 'tau'), 'lex')
ctxzw = flint.fmpz_mpoly_ctx.get(('z', 'w'), 'lex')
z, w = ctxzw.gens()
sub_nu, sub_tau = SNU + 2 * z + w, S3 + w
npos = nneg = nterms = 0
c0 = None
worst = []
minpos = None
for (a, b), dct in sorted(s2.items()):
    q = ctxnt.from_dict(dct).compose(sub_nu, sub_tau)
    for e, c in q.to_dict().items():
        c = int(c)
        nterms += 1
        if c > 0:
            npos += 1
        elif c < 0:
            nneg += 1
            if len(worst) < 12:
                worst.append(((a, b) + tuple(e), c))
        if a == 0 and b == 0 and tuple(e) == (0, 0):
            c0 = c
print(f'(2) box expansion: {nterms} nonzero terms, positive {npos}, negative {nneg}; constant term sign',
      (c0 > 0) - (c0 < 0) if c0 is not None else None, f'({time.time()-T0:.1f}s)', flush=True)
if nneg:
    print('first negative terms ((x, y, z, w) exponents, digits of |coefficient|):', [(e, len(str(-c))) for e, c in worst])
ok = c0 is not None and c0 != 0 and (nneg == 0 or npos == 0) and ((c0 > 0) == (npos > 0))
print('P ONE-SIGNED ON THE BOX (all coefficients of one sign, nonzero constant term): PASS' if ok else 'NOT ONE-SIGNED')
