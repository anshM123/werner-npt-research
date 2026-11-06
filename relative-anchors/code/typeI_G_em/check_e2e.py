"""End-to-end exact checks of whole type I blocks of (R_inf)_{k-1} at rational weights with equal U-weights, independent
of the closed forms of equal_m.py.

For each configuration (T = copies 0..m-1 with arbitrary rational amplitudes, U = u copies of one amplitude om):
 1. the block from the Lemma B builder (fastj1.FastJ1, which applies E^{-1} sector by sector; the closed form cf.py is
    not used) equals the closed form cf.TypeI (both in units of 2Y), entry by entry;
 2. Y (p (x) f_h) = 0 for every h in U (the kernel vectors);
 3. with the coordinates P_{r0,h} deleted (r0 the first copy of U), all LDL^T pivots of the block are positive, i.e.
    the block is >= 0 with kernel exactly span{p (x) f_h}; the Schur complement G on the letter-2 states is > 0;
 4. the five weighted block Gershgorin quantities (H = 1/2), computed exactly from this G (emlib.direct_conditions),
    are positive, and equal the values of the very rational functions certified in prove_em.py at the corresponding
    point (a, b, r, s, x) = (s_t, s_t', sigma - s_t - s_t', om^2, u - 3) (for B: (s_t, 0, sigma - s_t, om^2, u - 3)).
python check_e2e.py"""
import sys, os, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
import flint
from cf import TypeI
from emlib import direct_G, direct_conditions, ldl_pivots, qf, Q
import prove_em as P

H = Q(1, 2)
_, _, CONDS = P.run(H)


def val(rf, pt):
    pt = [Q(v) for v in pt]
    return rf.n(*pt) / rf.d(*pt)


def checks(wT, om, u):
    m = len(wT)
    w = [Fr(x) for x in wT] + [Fr(om)] * u
    k = len(w)
    T = list(range(m))
    U = list(range(m, k))
    Gz, piv, Mb, lab = direct_G(wT, om, u)
    Mc = TypeI(w, T).build()
    n = len(lab)
    same = all(Mb[i][j] == Mc[i][j] for i in range(n) for j in range(n))
    idx = {l: i for i, l in enumerate(lab)}
    ker = True
    for h in U:
        v = [Fr(0)] * n
        for r in range(k):
            v[idx[('P', r, h)]] = w[r]
        ker &= all(sum(Mb[i][j] * v[j] for j in range(n)) == 0 for i in range(n))
    pd = all(p > 0 for p in piv)
    Gpd = all(p > 0 for p in ldl_pivots(flint.fmpq_mat(Gz)))
    dc = direct_conditions(wT, om, u, H, Gz)
    pos = all(v > 0 for vals in dc.values() for v in vals)
    s = qf(Fr(om)) ** 2
    sT = [qf(Fr(x)) ** 2 for x in wT]
    sig = sum(sT, Q(0))
    match = True
    ia = 0
    for t in range(m):
        ptB = (sT[t], 0, sig - sT[t], s, u - 3)
        for name in ('B_inv_00', 'B_inv_det', 'B_std'):
            match &= val(CONDS[name], ptB) == dc[name][t]
        for tp in range(m):
            if tp == t:
                continue
            ptA = (sT[t], sT[tp], sig - sT[t] - sT[tp], s, u - 3)
            for name in ('A_inv', 'A_std'):
                match &= val(CONDS[name], ptA) == dc[name][ia]
            ia += 1
    return same, ker, pd, Gpd, pos, match


rng = random.Random(11052026)
amps = [Fr(1, 9), Fr(1, 4), Fr(1, 2), Fr(2, 3), Fr(1), Fr(3, 2), Fr(2), Fr(3), Fr(7, 2), Fr(1, 30), Fr(12)]
ok = True
t0 = time.time()
cnt = 0
for (u, m, reps) in ((3, 2, 4), (3, 3, 3), (3, 4, 2), (3, 5, 1), (4, 2, 3), (4, 3, 2), (4, 4, 1), (5, 2, 2), (5, 3, 1),
                     (6, 2, 1), (7, 2, 1)):
    for rep in range(reps):
        wT = [rng.choice(amps) for _ in range(m)]
        om = rng.choice(amps)
        res = checks(wT, om, u)
        good = all(res)
        ok &= good
        cnt += 1
        print('u=%d m=%d wT=%s om=%s: Lemma B = closed form %s, kernel %s, pivots > 0 %s, G > 0 %s, Gershgorin '
              'quantities > 0 %s, = certified functions %s [%.1fs]' % (u, m, [str(x) for x in wT], om, *res,
                                                                     time.time() - t0), flush=True)
print('%d configurations' % cnt)
print('PASS' if ok else 'FAIL')
