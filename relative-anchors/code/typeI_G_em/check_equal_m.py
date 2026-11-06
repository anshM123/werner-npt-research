"""Exact validation of the symbolic equal-U-weight blocks of equal_m.py against direct exact elimination of the whole
type I block (cf.TypeI over Q, gauge P_{r0,h} deleted), for u = 3..7 and m = 1..4.

For each rational configuration (T-amplitudes arbitrary, all U-amplitudes equal to om) the exact G_z (z-coordinates
zeta_{t,h} = w_t w_h z_{t,h}) is projected on the invariant and standard parts:
  inv block (t, t'):  [[G_z[rho_t, rho_t'], sum_h' G_z[rho_t, z_t'h']], [sum_h G_z[z_th, rho_t'], sum_{h,h'} G_z[z_th, z_t'h']]]
  std block (t, t'):  G_z[z_th, z_t'h] - G_z[z_th, z_t'h']  (h != h'),
and compared with the closed forms evaluated at (u, sigma, s, s_t, s_t') (diagonal blocks at s1 = s_t).  Also checks:
the diagonal closed forms do not involve s2 (resp. s1), every entry of G_z is equivariant (the exact G_z is
invariant under the permutations of U, so the projections are complete), and a second gauge (P_{r1,h} deleted, r1 the
second copy of U) gives the same G_z at some points."""
import sys, os, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from emlib import load, ev, Gz, Q

R = load()
ok = True
# the diagonal blocks depend only on their own T-weight
for key in ('inv', 'std'):
    G = R[key]
    nk = len(G) // 2
    for i in range(nk):
        for j in range(nk):
            for (n, d) in (G[i][j],):
                for p in (n, d):
                    if any(m[4] != 0 for m in p.monoms()):
                        ok = False
                        print('diagonal block (t1) depends on s2:', key, i, j)
            for (n, d) in (G[nk + i][nk + j],):
                for p in (n, d):
                    if any(m[3] != 0 for m in p.monoms()):
                        ok = False
                        print('diagonal block (t2) depends on s1:', key, i, j)
print('diagonal blocks depend only on their own T-weight:', ok)

rng = random.Random(20261105)
amps = [Q(1, 7), Q(1, 3), Q(1, 2), Q(2, 3), Q(1), Q(3, 2), Q(2), Q(5, 2), Q(4), Q(1, 20), Q(10)]
t0 = time.time()
npts = 0
nblocks = 0
for u in range(3, 8):
    for m in range(1, 5):
        if u + m > 9:
            continue
        reps = 3 if u + m <= 7 else 2
        for rep in range(reps):
            wT = [rng.choice(amps) for _ in range(m)]
            om = rng.choice(amps)
            wU = [om] * u
            Gq, mm, uu = Gz(wT, wU)
            d = u + 1
            sT = [x * x for x in wT]
            sig = sum(sT, Q(0))
            s = om * om
            good = True
            # equivariance (permutation (h0 h1) of U) of G_z
            perm = list(range(d))
            perm[1], perm[2] = perm[2], perm[1]
            P = [t * d + perm[i] for t in range(m) for i in range(d)]
            for i in range(m * d):
                for j in range(m * d):
                    if Gq[i, j] != Gq[P[i], P[j]]:
                        good = False
            for t in range(m):
                for tp in range(m):
                    if t == tp:
                        vals = [u, sig, s, sT[t], 0]
                        inv = [[ev(R['inv'][i][j], vals) for j in range(2)] for i in range(2)]
                        std = ev(R['std'][0][0], vals)
                    else:
                        vals = [u, sig, s, sT[t], sT[tp]]
                        inv = [[ev(R['inv'][i][2 + j], vals) for j in range(2)] for i in range(2)]
                        std = ev(R['std'][0][1], vals)
                    b0, b1 = t * d, tp * d
                    e00 = Gq[b0, b1]
                    e01 = sum((Gq[b0, b1 + 1 + h] for h in range(u)), Q(0))
                    e10 = sum((Gq[b0 + 1 + h, b1] for h in range(u)), Q(0))
                    e11 = sum((Gq[b0 + 1 + h, b1 + 1 + hp] for h in range(u) for hp in range(u)), Q(0))
                    es = Gq[b0 + 1, b1 + 1] - Gq[b0 + 1, b1 + 2]
                    if [e00, e01, e10, e11, es] != [inv[0][0], inv[0][1], inv[1][0], inv[1][1], std]:
                        good = False
                        print('  MISMATCH u=%d m=%d t=%d tp=%d' % (u, m, t, tp), flush=True)
                    nblocks += 1
            if rep == 0 and u + m <= 7:
                G2, _, _ = Gz(wT, wU, gauge_r0=m + 1)
                if G2 != Gq:
                    good = False
                    print('  gauge dependence u=%d m=%d' % (u, m))
            npts += 1
            ok &= good
            print('u=%d m=%d wT=%s om=%s: %s [%.1fs]' % (u, m, [str(x) for x in wT], om, 'equal' if good else 'MISMATCH',
                                                       time.time() - t0), flush=True)
print('%d configurations, %d (t, t\') blocks compared (invariant 2 x 2 and standard 1 x 1)' % (npts, nblocks))
print('PASS' if ok else 'FAIL')
