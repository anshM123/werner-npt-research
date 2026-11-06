"""Exploration (exact data, 120-digit linear algebra; no claim): the proof plan
    Psi > 0  <=  Psi_pp >= kappa D_pp  and  Psi_gg - Psi_gp (kappa D_pp)^{-1} Psi_gp^T > 0,
where g = (gamma, xi, rho) are the global directions (ch = gamma 1, zh = xi 1, rho) and p the perpendicular subspace
{sum s_h ah_h = sum s_h ch_h = sum s_h zh_h = 0, rho = 0}, D_pp the local part restricted to p.  Reports, for
kappa = 0.65 and 0.85, the smallest eigenvalue of the 3 x 3 matrix after the congruence by diag(Psi_gg)^{-1/2}."""
import sys, random
import mpmath as mp
import flint
from reduced_c import ReducedC, assemble
from certify_points import family

mp.mp.dps = 120
R = ReducedC()
rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 5)


def tomp(Mq):
    A = mp.matrix(Mq.nrows(), Mq.ncols())
    for i in range(Mq.nrows()):
        for j in range(Mq.ncols()):
            A[i, j] = mp.mpf(int(Mq[i, j].p)) / int(Mq[i, j].q)
    return A


for u in [int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else '3,4,6,8').split(',')]:
    worst = {0.65: 9.0, 0.85: 9.0}
    for k in range(8):
        wt, wU = family(u, k, rng)
        num = R.numeric(wt, wU)
        Psi, Dm, Um, Cm = assemble(num)
        n = 3 * u + 1
        sv = num['sv']
        gcols = []
        for blk in (1, 2):
            v = [flint.fmpq(0)] * n
            for h in range(u):
                v[blk * u + h] = flint.fmpq(1)
            gcols.append(v)
        v = [flint.fmpq(0)] * n
        v[n - 1] = flint.fmpq(1)
        gcols.append(v)
        pcols = []
        for blk in range(3):
            for h in range(1, u):
                v = [flint.fmpq(0)] * n
                v[blk * u + h] = flint.fmpq(1)
                v[blk * u] = -sv[h] / sv[0]
                pcols.append(v)
        Tg = flint.fmpq_mat(n, 3, [gcols[j][i] for i in range(n) for j in range(3)])
        m = len(pcols)
        Tp = flint.fmpq_mat(n, m, [pcols[j][i] for i in range(n) for j in range(m)])
        Pgg = tomp(Tg.transpose() * Psi * Tg)
        Pgp = tomp(Tg.transpose() * Psi * Tp)
        Dpp = tomp(Tp.transpose() * Dm * Tp)
        for kap in (0.65, 0.85):
            X = mp.inverse(Dpp * kap) * Pgp.T
            Mg = Pgg - Pgp * X
            sc = [1 / mp.sqrt(Pgg[i, i]) for i in range(3)]
            Mn = mp.matrix(3, 3)
            for i in range(3):
                for j in range(3):
                    Mn[i, j] = Mg[i, j] * sc[i] * sc[j]
            ev = mp.eigsy((Mn + Mn.T) / 2, eigvals_only=True)
            worst[kap] = min(worst[kap], min(float(e) for e in ev))
    print(f'u = {u}: smallest normalised eigenvalue of the global 3 x 3 condition: kappa = 0.65: {worst[0.65]:.3e}; '
          f'kappa = 0.85: {worst[0.85]:.3e}', flush=True)
