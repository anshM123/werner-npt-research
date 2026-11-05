"""Sign certificates for the level reductions of lev.py.

python cert_lev.py FORM l off_1 ... off_l [bscale] [order=i,j,...] [only=inv|std]
  FORM = F: the invariant form Finv and the standard forms Fstd of the P-first reduction Phi;
  FORM = G: the gauge-fixed block (state P_t deleted), Ginv and Gstd, with the B-part scaled by bscale.
  order: a permutation of the invariant coordinates (default: the stored order (rho_t, rho_1..rho_l, q_1..q_l) for F).
Level sizes u_i = off_i + x_i, x_i >= 0.  Every leading principal minor (in the chosen order) is factored over Q; it
passes if every irreducible factor of numerator and denominator has coefficients of one sign and a monomial free of the
x_i, and the resulting sign is positive.  Then the form is positive definite for all p, q_i > 0 and all real
u_i >= off_i (Sylvester).  The standard form of level i is only certified if off_i >= 2 (it exists only for u_i >= 2).
Stops at the first minor that does not pass.  Prints PASS if every minor passes."""
import sys, time
from fractions import Fraction as Fr
from lev import build, det, signcheck

t0 = time.time()
args = [a for a in sys.argv[1:] if '=' not in a]
opts = dict(a.split('=', 1) for a in sys.argv[1:] if '=' in a)
form = args[0]
l = int(args[1])
offs = [int(a) for a in args[2:2 + l]]
bsc = Fr(args[2 + l]) if len(args) > 2 + l else Fr(1)
fixed = [v == '1' for v in opts['fixed'].split(',')] if 'fixed' in opts else None
S = build(l, offs, bsc, want_phi=(form == 'F'), fixed=fixed)
K = S['K']
inv, std = (S['Finv'], S['Fstd']) if form == 'F' else (S['Ginv'], S['Gstd'])
if 'order' in opts:
    o = [int(v) for v in opts['order'].split(',')]
    assert sorted(o) == list(range(len(inv)))
    inv = [[inv[i][j] for j in o] for i in o]
print(f'{form}, l = {l}, offsets {offs}, fixed {fixed}, bscale {bsc}, order {opts.get("order", "stored")}: invariant form of size '
      f'{len(inv)} [{time.time() - t0:.1f}s]', flush=True)
ok = True
mats = []
if opts.get('only', 'inv') == 'inv' or 'only' not in opts:
    mats.append(('invariant', inv))
if opts.get('only', 'std') == 'std' or 'only' not in opts:
    mats += [(f'standard, level {i + 1}', std[i]) for i in range(l) if offs[i] >= 2]
for name, A in mats:
    sym = all(A[i][j] == A[j][i] for i in range(len(A)) for j in range(len(A)))
    ok &= sym
    for m in range(1, len(A) + 1):
        dm = det([row[:m] for row in A[:m]], K(0))
        sg, desc = signcheck(dm)
        good = sg == 1
        ok &= good
        print(f'  {name} (symmetric {sym}) minor {m}: sign {sg}; numerator deg {dm.n.total_degree()} '
              f'({len(dm.n)} terms); factors (deg/terms^mult): {desc} [{time.time() - t0:.1f}s]', flush=True)
        if not good:
            break
    if not ok:
        break
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
