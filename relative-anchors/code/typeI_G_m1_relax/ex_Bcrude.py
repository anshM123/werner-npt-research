"""Exploration (exact data; no claim): the crude version of (B_k) without the constraint of p:
  Psi_gg - (1/kappa) sum_h b_h^T D_h^{-1} b_h > 0,   b_h = the coupling of g to the local variables of copy h
(full-space representative), kappa = 1/2; versus the exact constrained version (B_k).  Positivity by exact minors."""
import sys, random
import flint
from reduced_c import ReducedC, assemble
from certify_points import family, leading_minors_bareiss
Q = flint.fmpq
R = ReducedC()
rng = random.Random(11)
def pd(Mx):
    mins, _ = leading_minors_bareiss(Mx)
    return len(mins) == Mx.nrows() and all(x > 0 for x in mins)
for u in (4, 6, 8):
    good = 0
    fails = []
    for k in range(8):
        wt, wU = family(u, k, rng)
        num = R.numeric(wt, wU)
        Psi, Dm, Um, Cm = assemble(num)
        n = 3 * u + 1
        gcols = []
        for blk in (1, 2):
            v = [Q(0)] * n
            for h in range(u):
                v[blk * u + h] = Q(1)
            gcols.append(v)
        v = [Q(0)] * n; v[n - 1] = Q(1); gcols.append(v)
        Tg = flint.fmpq_mat(n, 3, [gcols[j][i] for i in range(n) for j in range(3)])
        Pgg = Tg.transpose() * Psi * Tg
        Cg = Tg.transpose() * Psi          # 3 x n: coupling of g to every coordinate
        corr = flint.fmpq_mat(3, 3)
        for h in range(u):
            idx = [h, u + h, 2 * u + h]
            bh = flint.fmpq_mat(3, 3, [Cg[i, j] for i in range(3) for j in idx])
            Dh = flint.fmpq_mat(3, 3, [Dm[i, j] for i in idx for j in idx])
            corr += bh * Dh.solve(bh.transpose())
        ok = pd(Pgg - corr * 2)
        good += ok
        if not ok:
            fails.append(k % 8)
    print(f'u = {u}: crude (B_1/2) holds at {good}/8; fails at families {fails}', flush=True)
