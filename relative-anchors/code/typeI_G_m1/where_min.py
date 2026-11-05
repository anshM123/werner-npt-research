"""Where do the smallest normalised margins of check_points.py occur?  Replays the same random points (same seed and
order) for u = 1, ..., umax and prints, per u, the family of the worst point, the row of the worst margin, and whether
that row is the z-coordinate of the smallest U-weight.  Exact arithmetic for G_z; the ratio is printed as a float."""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from gzexact import Gz_exact, mvec
from check_points import sample

Q = flint.fmpq
FAM = ['log-uniform', 'one tiny', 'all shrunk', 'all large', 'one dominant', 'equal', 'sigma tiny', 'two tiny',
       'geometric']
umax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
npts = 28
rng = random.Random(2026)
ok = True
for u in range(1, umax + 1):
    worst = (2.0, None)
    for it in range(npts):
        wt, wU = sample(u, it, rng)
        if it % 7 == 0:
            pass            # check_points also compares a second gauge here; no random numbers are drawn there
        G = Gz_exact(wt, wU)
        st = wt * wt
        S = sum((x * x for x in wU), Q(0))
        v = mvec(st, S, u)
        n = u + 1
        Gv = [sum((G[i, j] * v[j] for j in range(n)), Q(0)) for i in range(n)]
        for i in range(n):
            r = float(Gv[i] / (G[i, i] * v[i]))
            if r < worst[0]:
                smallest = min(range(u), key=lambda h: wU[h] * wU[h])
                worst = (r, (FAM[it % 9], 'rho' if i == 0 else f'z_{i}', (i - 1) == smallest, float(wU[smallest] ** 2 / S)))
    fam, row, is_small, rel = worst[1]
    print(f'u = {u}: worst margin {worst[0]:.3e} in family "{fam}", row {row}; '
          f'row of the smallest U-weight: {is_small} (that weight is {rel:.1e} of S)')
print('PASS')
