"""Exact check of the per-copy form of (B_k) (NOTES.md, Section 13, Lemma P2): with b_h the 3 x 3 coupling of the global
directions g to the local variables of copy h (full-space representative), D_h the local block, and
  B0 = sum_h b_h^T D_h^{-1} b_h,  B1 = sum_h s_h D_h^{-1} b_h,  B2 = sum_h s_h^2 D_h^{-1},
the coupling term of (B_k) on the perpendicular subspace equals (1/kappa)(B0 - B1^T B2^{-1} B1):
  Psi_gp (kappa D_pp)^{-1} Psi_gp^T = (1/kappa)(B0 - B1^T B2^{-1} B1)   (exact equality at rational points)."""
import random
import flint
from reduced_c import ReducedC, assemble
from certify_points import family
Q = flint.fmpq
R = ReducedC()
rng = random.Random(41)
ok = True
for u in (3, 4, 6, 8):
    for k in range(4):
        wt, wU = family(u, k, rng)
        num = R.numeric(wt, wU)
        Psi, Dm, Um, Cm = assemble(num)
        n = 3 * u + 1
        sv = num['sv']
        gcols = []
        for blk in (1, 2):
            v = [Q(0)] * n
            for h in range(u):
                v[blk * u + h] = Q(1)
            gcols.append(v)
        v = [Q(0)] * n; v[n - 1] = Q(1); gcols.append(v)
        Tg = flint.fmpq_mat(n, 3, [gcols[j][i] for i in range(n) for j in range(3)])
        pcols = []
        for blk in range(3):
            for h in range(1, u):
                v = [Q(0)] * n
                v[blk * u + h] = Q(1)
                v[blk * u] = -sv[h] / sv[0]
                pcols.append(v)
        m = len(pcols)
        Tp = flint.fmpq_mat(n, m, [pcols[j][i] for i in range(n) for j in range(m)])
        Pgp = Tg.transpose() * Psi * Tp
        Dpp = Tp.transpose() * Dm * Tp
        lhs = Pgp * Dpp.solve(Pgp.transpose())
        Cg = Tg.transpose() * Psi
        B0 = flint.fmpq_mat(3, 3); B1 = flint.fmpq_mat(3, 3); B2 = flint.fmpq_mat(3, 3)
        for h in range(u):
            idx = [h, u + h, 2 * u + h]
            bh = flint.fmpq_mat(3, 3, [Cg[j, i] for i in idx for j in range(3)])    # rows: classes, cols: g
            Dh = flint.fmpq_mat(3, 3, [Dm[i, j] for i in idx for j in idx])
            Dib = Dh.solve(bh)
            B0 += bh.transpose() * Dib
            B1 += Dib * sv[h]
            B2 += Dh.inv() * (sv[h] * sv[h])
        rhs = B0 - B1.transpose() * B2.solve(B1)
        eq = all(lhs[i, j] == rhs[i, j] for i in range(3) for j in range(3))
        ok &= eq
    print(f'u = {u}: Psi_gp D_pp^-1 Psi_gp^T == B0 - B1^T B2^-1 B1 exactly at 4 weight vectors: {ok}', flush=True)
print('PASS' if ok else 'FAIL')
