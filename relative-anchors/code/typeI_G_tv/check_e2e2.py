"""End-to-end exact checks of whole type I blocks of (R_inf)_{k-1} at two-valued weight vectors (squared weights s and
1 on U and on T), independent of the closed forms of tl_m.py.

For each configuration (U = u1 copies of amplitude om + u2 copies of amplitude 1; T = m1 copies of amplitude om + m2
copies of amplitude 1):
 1. the block from the Lemma B builder (fastj1.FastJ1, E^{-1} sector by sector) equals the closed form cf.TypeI;
 2. Y (p (x) f_h) = 0 for every h in U;
 3. with the coordinates P_{r0,h} deleted, all LDL^T pivots are positive (the block is >= 0 with kernel exactly
    span{p (x) f_h}); the Schur complement G on the letter-2 states is > 0;
 4. the weighted block Gershgorin quantities of prove_gersh.py (H = 1/2), computed exactly from this G (invariant part
    in the frame (rho, y^1, y^2), standard parts of both levels), are all positive; for u1 = 1 they are compared with the
    very rational functions certified by prove_gersh.py at the same point (x2, r, s).
python check_e2e2.py"""
import sys, os, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
import flint
from cf import TypeI
from emlib import lemmab_block, ldl_pivots, qf

Q = flint.fmpq
H = Q(1, 2)


def direct(wT, wU, bscale=1):
    m, u = len(wT), len(wU)
    w = [Fr(x) for x in list(wT) + list(wU)]
    k = len(w)
    T, U = list(range(m)), list(range(m, k))
    Mb, lab = lemmab_block(w, T, bscale)
    n = len(lab)
    idx = {l: i for i, l in enumerate(lab)}
    r0 = U[0]
    keep = [i for i in range(n) if not (lab[i][0] == 'P' and lab[i][1] == r0)]
    piv = ldl_pivots(flint.fmpq_mat([[qf(Mb[i][j]) for j in keep] for i in keep]))
    tb = []
    for t in T:
        tb.append(idx[('R', t)])
        for h in U:
            tb.append(idx[('Z', t, h)])
    el = [i for i in keep if i not in tb]
    C = flint.fmpq_mat([[qf(Mb[i][j]) for j in el] for i in el])
    B = flint.fmpq_mat([[qf(Mb[i][j]) for j in tb] for i in el])
    A = flint.fmpq_mat([[qf(Mb[i][j]) for j in tb] for i in tb])
    G = A - B.transpose() * C.solve(B)
    d = []
    for t in T:
        d.append(Q(1))
        for h in U:
            d.append(qf(w[t]) * qf(w[h]))
    Gz = [[G[i, j] * d[i] * d[j] for j in range(len(d))] for i in range(len(d))]
    return Gz, piv, Mb, lab, w, T, U


def quantities(Gz, sT, s, u1, u2):
    """Gershgorin quantities (dict name -> list of values) from the exact G_z; level 1 = copies 1..u1 (weight s)."""
    m = len(sT)
    dd = 1 + u1 + u2
    sig = sum(sT, Q(0))
    S1, S2 = u1 * s, Q(u2)
    S = S1 + S2
    I = [[0], list(range(1, 1 + u1)), list(range(1 + u1, dd))]
    Nm = [[Q(1), Q(0), Q(0)], [Q(0), 1 / S, S2 / S], [Q(0), 1 / S, -S1 / S]]

    def inv_blk(t, tp):
        b0, b1 = t * dd, tp * dd
        Gz3 = [[sum((Gz[b0 + i][b1 + j] for i in I[p] for j in I[q]), Q(0)) for q in range(3)] for p in range(3)]
        return [[sum((Nm[k][i] * Gz3[k][l] * Nm[l][j] for k in range(3) for l in range(3)), Q(0)) for j in range(3)]
                for i in range(3)]

    def std_blk(t, tp, lvl):
        b0, b1 = t * dd, tp * dd
        h, hp = I[lvl][0], I[lvl][1]
        return Gz[b0 + h][b1 + h] - Gz[b0 + h][b1 + hp]

    out = {'A_inv': [], 'B_inv': [], 'A_std1': [], 'B_std1': [], 'A_std2': [], 'B_std2': []}
    kap = lambda st: [Q(1), st * sig / (S * (sig + S)), st * S1 * S2 / S]
    sl = {1: s, 2: Q(1)}
    for t in range(m):
        th1 = 1 - sT[t] / sig
        g = inv_blk(t, t)
        kt = kap(sT[t])
        M = [[g[i][j] - (H * th1 * kt[i] if i == j else 0) for j in range(3)] for i in range(3)]
        out['B_inv'].append(ldl_pivots(flint.fmpq_mat(M)))
        for lvl, uu in ((1, u1), (2, u2)):
            if uu >= 2:
                out['B_std%d' % lvl].append(std_blk(t, t, lvl) / (sT[t] * sl[lvl]) - H * th1)
        for tp in range(m):
            if tp == t:
                continue
            g = inv_blk(t, tp)
            ka, kb = kap(sT[t]), kap(sT[tp])
            fro = sum((g[i][j] ** 2 / (ka[i] * kb[j]) for i in range(3) for j in range(3)), Q(0))
            tgt = H * H * sT[t] * sT[tp] / sig ** 2
            out['A_inv'].append((t, tp, tgt - fro))
            for lvl, uu in ((1, u1), (2, u2)):
                if uu >= 2:
                    out['A_std%d' % lvl].append((t, tp, tgt - std_blk(t, tp, lvl) ** 2 / (sT[t] * sT[tp] * sl[lvl] ** 2)))
    return out


def main():
    import prove_gersh as PG
    t0 = time.time()
    conds = {name: (val, free) for name, val, free in PG.conditions('one', H)}
    print('certified functions for u1 = 1 built [%.1fs]' % (time.time() - t0), flush=True)
    rng = random.Random(6112026)
    amps = [Fr(1, 9), Fr(1, 4), Fr(1, 2), Fr(2, 3), Fr(3, 2), Fr(2), Fr(3), Fr(1, 30), Fr(12)]
    ok = True
    cnt = 0
    for (u1, u2, m1, m2, reps) in ((1, 2, 0, 2, 2), (1, 2, 1, 1, 2), (1, 2, 2, 0, 1), (1, 3, 1, 2, 2), (1, 3, 2, 1, 1),
                                   (1, 4, 0, 3, 1), (1, 4, 1, 1, 1), (1, 5, 2, 0, 1), (2, 2, 1, 1, 2), (2, 2, 0, 2, 1),
                                   (2, 3, 2, 1, 1), (3, 2, 1, 1, 1), (2, 2, 2, 1, 1), (3, 3, 1, 1, 1)):
        for rep in range(reps):
            om = rng.choice(amps)
            wT = [om] * m1 + [Fr(1)] * m2
            wU = [om] * u1 + [Fr(1)] * u2
            Gz, piv, Mb, lab, w, T, U = direct(wT, wU)
            Mc = TypeI(w, T).build()
            n = len(lab)
            same = all(Mb[i][j] == Mc[i][j] for i in range(n) for j in range(n))
            idx = {l: i for i, l in enumerate(lab)}
            ker = True
            for h in U:
                v = [Fr(0)] * n
                for r in range(len(w)):
                    v[idx[('P', r, h)]] = w[r]
                ker &= all(sum(Mb[i][j] * v[j] for j in range(n)) == 0 for i in range(n))
            pd = all(p > 0 for p in piv)
            Gpd = all(p > 0 for p in ldl_pivots(flint.fmpq_mat(Gz)))
            s = qf(om) ** 2
            sT = [qf(x) ** 2 for x in wT]
            qn = quantities(Gz, sT, s, u1, u2)
            pos = all(all(p > 0 for p in pv) for pv in qn['B_inv']) and all(v > 0 for v in qn['B_std1'] + qn['B_std2'])
            pos &= all(v > 0 for (_, _, v) in qn['A_inv'] + qn['A_std1'] + qn['A_std2'])
            match = True
            if u1 == 1:
                sig = sum(sT, Q(0))
                lab_w = lambda st: 's' if st == s and om != 1 else '1'
                for (t, tp, v) in qn['A_inv']:
                    ka, kb = lab_w(sT[t]), lab_w(sT[tp])
                    if (ka, kb) == ('1', 's'):
                        ka, kb = 's', '1'
                        t, tp = tp, t
                    f, _ = conds["tl_inv_one A_inv (s_t, s_t') = (%s, %s)" % (ka, kb)]
                    pt = [Q(0), Q(u2 - 2), sig - sT[t] - sT[tp], s]
                    match &= f.n(*pt) / f.d(*pt) == v
                for t in range(len(sT)):
                    ka = lab_w(sT[t])
                    pt = [Q(0), Q(u2 - 2), sig - sT[t], s]
                    for p in range(3):
                        f, _ = conds['tl_inv_one B_inv s_t = %s pivot %d' % (ka, p + 1)]
                        match &= f.n(*pt) / f.d(*pt) == qn['B_inv'][t][p]
            good = same and ker and pd and Gpd and pos and match
            ok &= good
            cnt += 1
            print('u1=%d u2=%d m1=%d m2=%d om=%s: Lemma B = closed form %s, kernel %s, pivots > 0 %s, G > 0 %s, '
                  'Gershgorin quantities > 0 %s, = certified functions %s [%.1fs]'
                  % (u1, u2, m1, m2, om, same, ker, pd, Gpd, pos, match if u1 == 1 else '(n/a)', time.time() - t0),
                  flush=True)
    print('%d configurations' % cnt)
    return ok


if __name__ == '__main__':
    print('PASS' if main() else 'FAIL')
