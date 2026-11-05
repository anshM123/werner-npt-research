"""Leading-minor sign certificates for family 2 (the |T| = 1 blocks of (R_inf)_2), part by part.

python cert_f2.py PART U0 U1 [tw=free|0|1] [bscale=r]
  PART in inv, std, two, anti, mix (irrep on level 0; mix on both levels);
  U0, U1: a fixed level size (integer) or 'o+' (every size u = o + x, x >= 0 real);
  tw: the squared weight of t: a free variable a (default), or that of level 0 (= 1) or level 1 (= q).
Squared weights 1 on level 0 and q on level 1.  Steps: build the reduced blocks over Q(a, q, x_i) (f2sym.py), check
that the Gram determinants of the test vectors are positive on the regime (so the test vectors are bases of the
multiplicity spaces), check F kappa = 0 for the kernel vectors kappa of the part, delete the gauge coordinates (the
test vectors on the diagonal states P_{c,c}), then
  1. the R-block (diagonal) has positive entries;
  2. after eliminating the R-states, the Schur complement on the P-vectors (gauge-fixed) followed by the Q-vectors
     has positive leading principal minors.
A minor (or entry) passes if its numerator and denominator are products of irreducible polynomials with one-signed
coefficients, each with a monomial free of the x_i, of total sign +1: then it is positive for a, q > 0, x_i >= 0.
So the gauge-fixed form is positive definite on the whole regime.  Prints PASS or FAIL."""
import sys, time
from fractions import Fraction as Fr
from flintfield import Field
from f2sym import Model2, reduce_part2
from l1sym import det, minv, mm, mT
from cert import signcheck
import t2lev

t0 = time.time()
args = [x for x in sys.argv[1:] if '=' not in x]
opts = dict(x.split('=', 1) for x in sys.argv[1:] if '=' in x)
part, s0, s1 = args[0], args[1], args[2]
tw = opts.get('tw', 'free')
bsc = Fr(opts.get('bscale', '1'))


def parse(spec):
    if spec.endswith('+'):
        return None, int(spec[:-1])
    return int(spec), None


(f0, o0), (f1, o1) = parse(s0), parse(s1)
names = (['a'] if tw == 'free' else []) + (['q'] if (f1 is None or f1 > 0 or tw == '1') else []) + \
    (['x0'] if f0 is None else []) + (['x1'] if f1 is None else [])
K = Field(names)
g = dict(zip(names, K.gens))
one = K(1)
q = g.get('q', K(Fr(5, 2).numerator) / K(Fr(5, 2).denominator))
U = [g['x0'] + o0 if f0 is None else K(f0), g['x1'] + o1 if f1 is None else K(f1)]
a = g['a'] if tw == 'free' else (one if tw == '0' else q)
M = Model2(one, a, [one, q], U, bsc)
red = reduce_part2(M, part, f0 if f0 is not None else o0, f1 if f1 is not None else o1)
print(f'family 2, part {part}, u0 = {s0}, u1 = {s1}, t-weight {tw}, bscale {bsc}, variables {names}: vectors '
      f'{red["names"]} [{time.time() - t0:.1f}s]', flush=True)
ok = True
for gname in ('GR', 'GP', 'GQ'):
    if gname in red and len(red[gname]):
        sg, desc = signcheck(det(red[gname], one))
        print(f'  det of the Gram matrix {gname}: sign {sg} ({desc})', flush=True)
        ok &= sg == 1
from f2sym import assemble
F, iR, iP, iQ = assemble(red)
n = len(F)
ok &= all(F[i][j] == F[j][i] for i in range(n) for j in range(n))
if bsc == 1:
    for kv in red['ker']:
        vec = [K(0)] * n
        for j, c in enumerate(kv):
            vec[iP[j]] = K(c)
        z = all(sum((F[i][j] * vec[j] for j in range(n)), K(0)) == 0 for i in range(n))
        print(f'  kernel vector {kv}: F kappa = 0 identically: {z}', flush=True)
        ok &= z
iPg = [iP[j] for j in range(len(iP)) if j not in red['gauge']]


def check(name, val):
    global ok
    sg, desc = signcheck(val)
    good = sg == 1
    ok &= good
    print(f'  {name}: sign {sg}; numerator deg {val.n.total_degree()} ({len(val.n)} terms), denominator deg '
          f'{val.d.total_degree()} ({len(val.d)} terms); {desc} [{time.time() - t0:.1f}s]', flush=True)
    return good


# 1. R-block
RR = [[F[i][j] for j in iR] for i in iR]
diag = all(RR[i][j] == 0 for i in range(len(iR)) for j in range(len(iR)) if i != j)
print(f'  R-block diagonal: {diag}', flush=True)
ok &= diag
for i in range(len(iR)):
    check(f'R-entry {i + 1}', RR[i][i])
# eliminate R (diagonal), then the leading minors of the rest in the order P (gauge-fixed), Q
rest = iPg + iQ
S1 = [[F[i][j] - sum((F[i][r] * F[r][j] / F[r][r] for r in iR if F[i][r] != 0 and F[r][j] != 0), K(0))
       for j in rest] for i in rest]
res = []
lab = ['P'] * len(iPg) + ['Q'] * len(iQ)
t2lev.leading_minors(S1, one, lambda m, d: res.append(check(f'minor {m} ({lab[m - 1]})', d)) or res[-1])
ok &= len(res) == len(rest) and all(res)
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
