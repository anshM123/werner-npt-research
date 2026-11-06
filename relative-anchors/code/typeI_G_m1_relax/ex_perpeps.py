"""Exploration (exact data, 120-digit eigenvalues; no claim): on the perpendicular subspace
{sum s_h ah_h = sum s_h ch_h = sum s_h zh_h = 0, rho = 0}, the largest generalized eigenvalue
eps* = max (-y^T C y)/(sum_h v_h^T D_h v_h).  eps* < 1 means the global part consumes less than the local part there."""
import sys, random
import mpmath as mp
import flint
from reduced_c import ReducedC, assemble
from certify_points import family

mp.mp.dps = 120
R = ReducedC()
rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 5)
for u in [int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else '4,6,8').split(',')]:
    worst = 0
    for k in range(8):
        wt, wU = family(u, k, rng)
        num = R.numeric(wt, wU)
        Psi, Dm, Um, Cm = assemble(num)
        n = 3 * u + 1
        # constraint basis: vectors v in R^{3u} (rho excluded) with sum s_h v_{x,h} = 0 for each class x
        sv = num['sv']
        cols = []
        for blk in range(3):
            for h in range(1, u):
                v = [flint.fmpq(0)] * n
                v[blk * u + h] = flint.fmpq(1)
                v[blk * u] = -sv[h] / sv[0]
                cols.append(v)
        m = len(cols)
        T = flint.fmpq_mat(n, m, [cols[j][i] for i in range(n) for j in range(m)])
        G = T.transpose() * (Um * Cm * Um.transpose()) * T
        Dl = T.transpose() * Dm * T
        Gm = mp.matrix(m, m)
        Dmm = mp.matrix(m, m)
        for i in range(m):
            for j in range(m):
                Gm[i, j] = mp.mpf(int(G[i, j].p)) / int(G[i, j].q)
                Dmm[i, j] = mp.mpf(int(Dl[i, j].p)) / int(Dl[i, j].q)
        # generalized eigenvalues of (-G, Dl): Cholesky of Dl
        L = mp.cholesky(Dmm)
        Li = mp.inverse(L)
        Mx = -(Li * Gm * Li.T)
        ev = mp.eigsy((Mx + Mx.T) / 2, eigvals_only=True)
        epsk = max(float(e) for e in ev)
        worst = max(worst, epsk)
    print(f'u = {u}: max over 8 weight vectors of eps* = {worst:.4f}', flush=True)
