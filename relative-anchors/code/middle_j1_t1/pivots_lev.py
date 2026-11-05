"""Pivot certificate for a level family, by symmetric Gaussian elimination over the rational-function field.

python pivots_lev.py l off fix which
  e.g. 3 1,1,2 1,1,0 inv: two single U-copies of arbitrary weight and u_3 = 2 + x_3 equal ones.
The matrix is the invariant part (which = inv) or the standard part of level N (which = stdN) of D Phi D (Lemma P),
built by lev.py with p = 1 (every entry is homogeneous of degree 0 in the squared weights); Phi itself is checked
against j1lib by check_lev.py.  Every pivot (ratio of consecutive leading principal minors) must be a ratio of
polynomials with coefficients of one sign, positive, each with a monomial free of the x-variables (lev.signcheck).
Then every leading principal minor is positive for all q_i > 0 and all real level sizes u_i >= off_i, and the matrix is
positive definite (Sylvester).  Prints PASS if all pivots pass."""
import sys, time
from lev import build, signcheck

t0 = time.time()
l = int(sys.argv[1])
offs = [int(v) for v in sys.argv[2].split(',')]
fix = [v == '1' for v in sys.argv[3].split(',')]
which = sys.argv[4]
S = build(l, offs, 1, want_phi=True, fixed=fix, p_one=True)
A = S['Finv'] if which == 'inv' else S['Fstd'][int(which[3:]) - 1]
n = len(A)
print(f'levels {offs}, fixed {fix}, {which}: matrix of size {n} over Q({",".join(S["names"])}) (p = 1) '
      f'[{time.time() - t0:.1f}s]', flush=True)
ok = all(A[i][j] == A[j][i] for i in range(n) for j in range(n))
Mx = [row[:] for row in A]
for c in range(n):
    pv = Mx[c][c]
    sg, desc = signcheck(pv)
    good = sg == 1
    ok &= good
    print(f'  pivot {c + 1}/{n}: numerator deg {pv.n.total_degree()} ({len(pv.n)} terms), denominator deg '
          f'{pv.d.total_degree()} ({len(pv.d)} terms); sign {sg} {desc[:90]} [{time.time() - t0:.1f}s]', flush=True)
    if not good:
        break
    for i in range(c + 1, n):
        if Mx[i][c] != 0:
            f = Mx[i][c] / pv
            for j in range(c + 1, n):
                if Mx[c][j] != 0:
                    Mx[i][j] = Mx[i][j] - f * Mx[c][j]
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
