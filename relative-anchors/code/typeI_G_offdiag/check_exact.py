"""Exact check (python-flint) of the closed form rawcf.raw_pair against R_raw of exact_offdiag.py
(exact_offdiag.raw_split: the whole type I block from cf.TypeI, then psi, q and p eliminated exactly).

Free amplitudes are rational points of the unit sphere (inverse stereographic projection), so S = 1 exactly and every
w_h is rational; letter-2 amplitudes are arbitrary positive rationals.  Every off-diagonal block (t, t') is compared
entry by entry in the coordinates (rho_t, zeta_{t,h}).  Includes equal letter-2 weights and spectator copies.
Usage: python check_exact.py seed npoints"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from exact_offdiag import raw_split
from rawcf import raw_pair

Q = flint.fmpq


def sphere_point(rng, u):
    if u == 1:
        return [Q(1)]
    while True:
        v = [Q(rng.randint(1, 30), rng.randint(31, 90)) for _ in range(u - 1)]
        n2 = sum((x * x for x in v), Q(0))
        if n2 < 1:
            break
    return [2 * x / (1 + n2) for x in v] + [(1 - n2) / (1 + n2)]


if __name__ == '__main__':
    seed, npts = int(sys.argv[1]), int(sys.argv[2])
    rng = random.Random(seed)
    nb = 0
    for it in range(npts):
        m = rng.randint(2, 4); u = rng.randint(1, 5)
        wU = sphere_point(rng, u)
        assert sum((x * x for x in wU), Q(0)) == 1
        scale = Q(2) ** rng.randint(-4, 4)
        wT = [scale * Q(rng.randint(1, 40), rng.randint(1, 40)) for _ in range(m)]
        if it % 3 == 0:
            wT[1] = wT[0]          # equal letter-2 weights
        G, R, info = raw_split(wT, wU)
        d = u + 1
        sU = [x * x for x in wU]
        sig = sum((x * x for x in wT), Q(0))
        for t in range(m):
            for tp in range(m):
                if t == tp:
                    continue
                M, st = raw_pair(sU, sU, wT[t] ** 2, wT[tp] ** 2, sig, wT[t], wT[tp], F=Q)
                for i in range(d):
                    for j in range(d):
                        ci = 1 if i == 0 else 1 / wU[i - 1]
                        cj = 1 if j == 0 else 1 / wU[j - 1]
                        assert M[i][j] * ci * cj == R[(t, tp)][i, j], ('mismatch', it, t, tp, i, j)
                nb += 1
        print(f'point {it}: m={m} u={u} wT={[str(x) for x in wT]} wU={[str(x) for x in wU]}: all off-diagonal '
              f'blocks equal', flush=True)
    print(f'{nb} off-diagonal blocks at {npts} rational points: closed form = R_raw exactly: PASS')
