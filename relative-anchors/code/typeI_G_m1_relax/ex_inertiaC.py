"""Exploration (exact data, eigenvalues in 300-digit arithmetic; no claim rests on this): inertia of the 31 x 31
matrix C of the reduced system, and of C restricted to the functionals of each variable class, at rational weights."""
import sys, random
import mpmath as mp
import flint
from reduced_c import ReducedC
from certify_points import family

mp.mp.dps = 300
R = ReducedC()
rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 3)
names = R.gens
for u in [int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else '4,8,12').split(',')]:
    for k in range(4):
        wt, wU = family(u, k, rng)
        num = R.numeric(wt, wU)
        J = len(num['C'])
        Cm = mp.matrix(J, J)
        for i in range(J):
            for j in range(J):
                x = num['C'][i][j]
                Cm[i, j] = mp.mpf(int(x.p)) / int(x.q)
        ev = mp.eigsy(Cm, eigvals_only=True)
        evs = sorted([float(e) for e in ev])
        mx = max(abs(e) for e in evs)
        neg = sum(1 for e in evs if e < -1e-40 * mx)
        zer = sum(1 for e in evs if abs(e) <= 1e-40 * mx)
        print(f'u = {u}, family {k}: C has {neg} negative, {zer} zero, {J - neg - zer} positive eigenvalues; '
              f'most negative {evs[0]:.3e}, largest {evs[-1]:.3e}', flush=True)
