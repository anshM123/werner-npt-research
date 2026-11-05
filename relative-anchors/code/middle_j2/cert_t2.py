"""Leading-minor sign certificates for family 1 (the |T| = 2 blocks of (R_inf)_2), through Lemma P2.

python cert_t2.py SPEC0 [SPEC1] [part=inv|std0|std1|all] [form=Phi|G] [bscale=r] [order=i,j,...] [tw=A,B]
  SPEC: a level size: an integer (fixed size) or 'o+' (every size u = o + x, x >= 0 real).  One SPEC: one U-level;
  two SPECs: two U-levels with squared weights 1 and q.  Squared weights a on t and b on tau (all free, > 0), or with
  tw=A,B (each 'free' or a level number 0/1) the squared weight of that level (all weights then take two values).
Builds the invariant form (orbits t, U_0[, U_1] of the R-states, tau, U_0[, U_1] of the R'-states, U_0[, U_1] of the
Q-states) and the standard forms (3 x 3 on (R, R', Q) for each level of size >= 2) of Phi = C - B^T diag(1/N) B over
Q(a, b[, q], x_i), and checks every leading principal minor: it passes if numerator and denominator are products of
irreducible polynomials with one-signed coefficients, each with a monomial free of the x_i, with total sign +1.  Then
the minor is positive for all a, b, q > 0 and x_i >= 0, and by Sylvester's criterion the form is positive definite on
the whole regime.  form=G: the gauge-fixed block (P_t deleted) instead (used with bscale for the negative control).
Prints PASS if every minor passes, FAIL otherwise."""
import sys, time
from fractions import Fraction as Fr
import t2lev
from t2lev import det
from cert import signcheck

t0 = time.time()
args = [x for x in sys.argv[1:] if '=' not in x]
opts = dict(x.split('=', 1) for x in sys.argv[1:] if '=' in x)
specs = args
l = len(specs)
fixed, offsets = [], []
for sp in specs:
    if sp.endswith('+'):
        fixed.append(None)
        offsets.append(int(sp[:-1]))
    else:
        fixed.append(int(sp))
        offsets.append(0)
part = opts.get('part', 'all')
form = opts.get('form', 'Phi')
bsc = Fr(opts.get('bscale', '1'))
tw = tuple(None if v == 'free' else int(v) for v in opts.get('tw', 'free,free').split(','))
alg, K, names = t2lev.symbolic(l, fixed=fixed, offsets=offsets, bscale=bsc, tw=tw)
one = K(1)
print(f'family 1, levels {specs}, T-weights {tw}, form {form}, bscale {bsc}, variables {names} [{time.time() - t0:.1f}s]', flush=True)
Finv, Fstd = t2lev.forms(alg, form)
print(f'  forms built: invariant {len(Finv)} x {len(Finv)} [{time.time() - t0:.1f}s]', flush=True)
todo = []
if part in ('all', 'inv'):
    todo.append(('inv', Finv))
for lv in range(l):
    sz = fixed[lv] if fixed[lv] is not None else offsets[lv]
    if sz >= 2 and part in ('all', f'std{lv}'):
        todo.append((f'std{lv}', Fstd[lv]))
ok = True
for name, F in todo:
    n = len(F)
    sym = all(F[i][j] == F[j][i] for i in range(n) for j in range(n))
    print(f'  {name}: size {n}, symmetric: {sym}', flush=True)
    ok &= sym
    order = list(range(n))
    if 'order' in opts and name == 'inv':
        order = [int(v) for v in opts['order'].split(',')]
    Fo = [[F[i][j] for j in order] for i in order]
    res = []

    def cb(m, dm):
        sg, desc = signcheck(dm)
        res.append(sg == 1)
        print(f'  {name} minor {m}: sign {sg}; numerator deg {dm.n.total_degree()} ({len(dm.n)} terms), denominator '
              f'deg {dm.d.total_degree()} ({len(dm.d)} terms); {desc} [{time.time() - t0:.1f}s]', flush=True)
        return sg == 1
    t2lev.leading_minors(Fo, one, cb)
    ok &= len(res) == n and all(res)
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
