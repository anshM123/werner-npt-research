"""Control for (B_k) of NOTES.md Section 13: with kappa much smaller than the true margin of (A), the coupling
correction (kappa D_pp)^{-1} is inflated; report at how many weight vectors (B_k) then fails (exact minors).
kappa = 1/20, 1/100, 1/1000."""
import random
import flint
from reduced_c import ReducedC, assemble
from certify_points import family, leading_minors_bareiss
Q = flint.fmpq
R = ReducedC()
def pd(Mx):
    mins, _ = leading_minors_bareiss(Mx)
    return len(mins) == Mx.nrows() and all(x > 0 for x in mins)
rng = random.Random(11)
for u in (4, 6):
    res = {Q(1, 20): 0, Q(1, 100): 0, Q(1, 1000): 0}
    for k in range(8):
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
        Pgg = Tg.transpose() * Psi * Tg
        Pgp = Tg.transpose() * Psi * Tp
        Dpp = Tp.transpose() * Dm * Tp
        X = Dpp.solve(Pgp.transpose())
        corr = Pgp * X
        for kap in res:
            res[kap] += (not pd(Pgg - corr * (1 / kap)))
    print(f'u = {u}: (B_kappa) fails at ' + ', '.join(f'{v}/8 for kappa = {k}' for k, v in res.items()), flush=True)
