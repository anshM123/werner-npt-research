"""Theorem G2 (the two-row part of the generic family [q0 >= 4, t >= 4, i >= 7], one marked copy) by the
determinant route: the assembly of the certificate.  It checks, exactly:
 1. the family form and its preliminaries (g2lib.build with check=True: N0 Q equals the numeric reduced form of
    mlib.py at 54 sample points, every entry denominator one-signed after the shifts, the split block-diagonal, the
    kernel vectors kappa_0, kappa_1 annihilated identically, the gauge determinant one-signed);
 2. the 12 irreducible factors f of the entry denominators of the gauge-fixed form Q (10 x 10) and, for each, the
    exponent e_f = max over permutations p of sum_a mult_f(den Q_{a,p(a)}) of D*; the candidate exponents of Den_c
    (g2slices.DEN_EXP) are at most e_f, so R = D*/Den_c is a polynomial;
 3. the degree bounds: T = D* det Q is a polynomial with deg_v T <= B_v = deg_v D* + max_p sum_a (deg_v num - deg_v den);
    the grid of the run is y = 0..B_t, z = 0..B_i and the low degrees are B_t - deg_t R, B_i - deg_i R;
 4. every factor of Den_c is one-signed (positive) after the shifts, so Den_c > 0 on the region;
 5. Q is positive definite at the point (q0, t, i, s) = (4, 4, 7, 1) (exact symmetric elimination);
 6. the run: every slice file of the grid exists, every column passed stage A (no coefficient above z^ZLOW), and
    stage B passed (no coefficient above y^YLOW, every coefficient of N~ nonnegative, a positive one free of x, y, z).
Conclusion (see NOTES): T and N~ R agree on the whole grid and have degrees <= B_t, B_i in y, z, so T = N~ R and
det Q = N~/Den_c > 0 on {x, y, z >= 0, s > 0}; with 5 and continuity, Q is positive definite on the whole region.
python g2final.py OUT"""
import sys, os, time, pickle
import flint
import g2core, g2lib, g2slices
from prove_fam import one_signed

out = sys.argv[1]
t0 = time.time()
ok = True
S = g2core.Slicer(check=True)
C, keep, M = S.C, S.keep, S.M
print(f'1. form and preliminaries checked [{time.time()-t0:.1f}s]', flush=True)

# 2. D* exponents
mult = {}
for a in keep:
    for b in keep:
        if M[a, b] == 0:
            continue
        c, fl = M[a, b].d.factor()
        for f, m in fl:
            mult[(str(f), a, b)] = m
estar = []
for f in S.factors:
    W = [[(mult.get((str(f), a, b), 0) if M[a, b] != 0 else None) for b in keep] for a in keep]
    estar.append(g2lib.maxperm(W))
dexp = g2slices.DEN_EXP
good2 = all(d <= e for d, e in zip(dexp, estar))
ok &= good2
print(f'2. factor exponents of D* {estar}, of Den_c {dexp}: Den_c divides D*: {good2}', flush=True)

# 3. degree bounds
deg = lambda v, exps: sum(e * f.degrees()[g2core.VAR[v]] for f, e in zip(S.factors, exps))
B = {}
for v in ('q0', 't', 'i', 's'):
    iv = g2core.VAR[v]
    W = [[(M[a, b].n.degrees()[iv] - M[a, b].d.degrees()[iv]) if M[a, b] != 0 else None for b in keep] for a in keep]
    B[v] = deg(v, estar) + g2lib.maxperm(W)
degR = {v: deg(v, estar) - deg(v, dexp) for v in ('q0', 't', 'i', 's')}
YMAX, ZMAX = B['t'], B['i']
YLOW, ZLOW = YMAX - degR['t'], ZMAX - degR['i']
print(f'3. degree bounds of T: {B}; degrees of R: {degR}; grid y = 0..{YMAX}, z = 0..{ZMAX}; '
      f'low degrees y^{YLOW}, z^{ZLOW}', flush=True)

# 4. Den_c > 0 on the region
signs = [one_signed(f, C) for f, e in zip(S.factors, dexp) if e > 0]
good4 = all(sg == 1 for sg in signs)
ok &= good4
print(f'4. every factor of Den_c positive after the shifts: {good4}', flush=True)

# 5. positive definite at one point
pt = [flint.fmpq(4), flint.fmpq(4), flint.fmpq(7), flint.fmpq(1)]
A = [[(M[a, b].n(*pt) / M[a, b].d(*pt)) if M[a, b] != 0 else flint.fmpq(0) for b in keep] for a in keep]
n = len(A)
pd = True
for c in range(n):
    p = A[c][c]
    if p <= 0:
        pd = False
        break
    for r in range(c + 1, n):
        f = A[r][c] / p
        for cc in range(c + 1, n):
            A[r][cc] -= f * A[c][cc]
ok &= pd
print(f'5. Q positive definite at (q0, t, i, s) = (4, 4, 7, 1): {pd}', flush=True)

# 6. the run
missing = [(y, z) for y in range(YMAX + 1) for z in range(ZMAX + 1)
           if not os.path.exists(os.path.join(out, f's_{y}_{z}.pkl'))]
colbad = []
for y in range(YMAX + 1):
    fn = os.path.join(out, f'col_{y}.pkl')
    if not os.path.exists(fn):
        colbad.append(y)
        continue
    with open(fn, 'rb') as fh:
        c = pickle.load(fh)
    if not (c['ok'] and c['zmax'] == ZMAX and c['zlow'] == ZLOW):
        colbad.append(y)
logB = open(os.path.join('logs', 'g2_stageB.log')).read().splitlines() if os.path.exists('logs/g2_stageB.log') else []
lines = [l for l in logB if l.startswith('N~(x, y, z, s)')]
stageb = bool(lines) and logB[-1].strip() == 'PASS' and ' 0 negative' in lines[-1] and \
    f'nonzero above y^{YLOW}: 0' in lines[-1] and 'free of x, y, z: True' in lines[-1]
good6 = not missing and not colbad and stageb
ok &= good6
print(f'6. slices missing: {len(missing)}; columns not passing stage A: {colbad}; stage B passed: {stageb}', flush=True)
if lines:
    print('   ', lines[-1])
print(f'total {time.time()-t0:.0f}s')
print('PASS' if ok else 'FAIL')
