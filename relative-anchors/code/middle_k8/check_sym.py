"""Checks of the certificates of cert_sym.py (isotypic reduction), exact, at rational squares q0:
 - every saved part F = B N B^T (N = den Yr) satisfies F(q0) = den(q0) B Yr(q0) B^T with Yr(q0) the rescaled block
   built by j1lib.py (Lemma B) at q0;
 - every saved leading minor of F equals the determinant of the leading principal submatrix of
   den(q0) B Yr(q0) B^T;
 - every basis vector of every part is orthogonal to the column-group indicators (and, for the letter-1 block of
   (R_inf)_{k-2}, to the indicator of the vacuum states), i.e. lies in W;
 - the parts are positive definite at q0 (exact symmetric elimination).
python check_sym.py FAMILIES [k=6]   (uses logs/cert_k<k>_<family>_<a>_<b>.pkl for the representative configurations)
Prints PASS or FAIL on the last line."""
import sys, time
from fractions import Fraction as Fr
import k5lib as K
from saved import load, ev_poly, ev_rf, cert_path

t0 = time.time()
args = [a for a in sys.argv[1:] if '=' not in a]
opts = dict(a.split('=', 1) for a in sys.argv[1:] if '=' in a)
k = int(opts.get('k', 6))
fams = args[0].split(',')
ok = True
for f in fams:
    j, m = K.fam(f)
    T = tuple(range(m))
    cnt = cm = 0
    for (a, b) in K.rep_configs(k, m):
        data = load(cert_path(f, a, b, k))
        st = data['states']
        n = len(st)
        groups = K.kernel_groups(k, st)
        vacs = [i for i, (a_, y) in enumerate(st) if not any(a_)] if (m == 0 and j == k - 2) else []
        for q0 in (Fr(4, 9), Fr(49, 16)):
            s = K.config_s(k, m, a, b, q0, Fr(1))
            w = [K.rsqrt(x) for x in s]
            st2, Mj = K.build_j1(k, j, T, w)
            assert st2 == st
            dv = ev_poly(data['den'], q0)
            for part in data['parts']:
                Bm = [[Fr(c) for c in row] for row in part['B']]
                d = len(Bm)
                inW = all(sum(row[i] for i in g) == 0 for row in Bm for g in groups.values()) and \
                    all(sum(row[i] for i in vacs) == 0 for row in Bm)
                MB = [[sum(Mj[i][l] * row[l] for l in range(n) if row[l] != 0) for i in range(n)] for row in Bm]
                R = [[dv * sum(Bm[r1][i] * MB[r2][i] for i in range(n) if Bm[r1][i] != 0) for r2 in range(d)]
                     for r1 in range(d)]
                F = [[ev_poly(x, q0) for x in row] for row in part['F']]
                eq = F == R
                mins = [ev_rf(x, q0) for x in part['minors']]
                meq = len(mins) == d and all(
                    Fr(int(K.to_fmpq_mat([row[:r] for row in R[:r]]).det().p),
                       int(K.to_fmpq_mat([row[:r] for row in R[:r]]).det().q)) == mins[r - 1] for r in range(1, d + 1))
                pd, _ = K.ldl_pd(R)
                good = inW and eq and meq and pd
                ok &= good
                cnt += 1
                cm += d
                if not good:
                    print(f'  {f} {(a, b)} q = {q0} part {part["lams"]}: in W {inW}, F = den B Y B^T {eq}, '
                          f'minors {meq}, PD {pd}')
    print(f'{f}: {cnt} (part, q) pairs, {cm} minors: parts equal the projections of the Lemma B block, minors equal '
          f'direct minors, parts PD: {ok} [{time.time() - t0:.1f}s]', flush=True)
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
