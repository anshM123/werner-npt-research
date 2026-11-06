"""Exact pointwise check of the two-part plan for Conjecture R (NOTES.md, Section 13), at rational weight vectors:
  (A_k)  Psi_pp - kappa D_pp >= 0 on the perpendicular subspace p (D_pp = the local part restricted to p);
  (B_k)  Psi_gg - Psi_gp (kappa D_pp)^{-1} Psi_gp^T > 0 on the global directions g = (gamma, xi, rho).
(A_k) and (B_k) together imply Psi > 0 (Schur complement, since Psi_pp^{-1} <= (kappa D_pp)^{-1}).  kappa = 1/2.
Also reports whether (A_k) fails for kappa = 99/100 (a control: the local part alone is not enough on p)."""
import sys, random
import flint
from reduced_c import ReducedC, assemble
from certify_points import family, leading_minors_bareiss

Q = flint.fmpq
R = ReducedC()
rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 11)
us = [int(x) for x in sys.argv[1].split(',')]
npts = int(sys.argv[2])
allok = True


def pd(Mx):
    mins, _ = leading_minors_bareiss(Mx)
    return len(mins) == Mx.nrows() and all(x > 0 for x in mins)


for u in us:
    cnt = dict(A=0, B=0, ctl=0)
    for k in range(npts):
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
        v = [Q(0)] * n
        v[n - 1] = Q(1)
        gcols.append(v)
        pcols = []
        for blk in range(3):
            for h in range(1, u):
                v = [Q(0)] * n
                v[blk * u + h] = Q(1)
                v[blk * u] = -sv[h] / sv[0]
                pcols.append(v)
        Tg = flint.fmpq_mat(n, 3, [gcols[j][i] for i in range(n) for j in range(3)])
        m = len(pcols)
        Tp = flint.fmpq_mat(n, m, [pcols[j][i] for i in range(n) for j in range(m)])
        Ppp = Tp.transpose() * Psi * Tp
        Dpp = Tp.transpose() * Dm * Tp
        Pgg = Tg.transpose() * Psi * Tg
        Pgp = Tg.transpose() * Psi * Tp
        kap = Q(1, 2)
        okA = pd(Ppp - Dpp * kap)
        X = (Dpp * kap).solve(Pgp.transpose())
        okB = pd(Pgg - Pgp * X)
        ctl = not pd(Ppp - Dpp * Q(99, 100))
        cnt['A'] += okA
        cnt['B'] += okB
        cnt['ctl'] += ctl
        allok &= okA and okB
    print(f'u = {u}: (A_1/2) holds at {cnt["A"]}/{npts}, (B_1/2) at {cnt["B"]}/{npts}; (A_99/100) fails at '
          f'{cnt["ctl"]}/{npts}', flush=True)
print('PASS' if allok else 'FAIL')
