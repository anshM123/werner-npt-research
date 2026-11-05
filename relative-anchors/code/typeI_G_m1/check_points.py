"""Exact rational-point checks of Conjecture M (m = 1): at rational weights w (s = w^2), G_z is computed exactly
(gzexact.py, gauge x_t = 0) and the following are checked in exact arithmetic:
  (Z) every off-diagonal entry of G_z is < 0;
  (D) G_z v > 0 entrywise, v = (4 s_t (s_t + S), 3 (2 s_t + S), ..., 3 (2 s_t + S));
  (P) the LDL^T pivots of G_z are all positive (G > 0 at the point; implied by (Z) and (D)).
At every 7th point it also checks that the gauge P_{r0,h} deleted gives the same G_z (gauge independence).
Families (cycled): log-uniform weights (spread up to e^{+-6}), one U-weight tiny (down to 1e-8 in s), all U-weights
shrunk (S/s_t down to 1e-8), all U-weights large (S/s_t up to 1e8), one U-weight dominant, equal U-weights, s_t tiny,
two U-weights tiny, geometric spread of the U-weights (ratio up to e per copy).
These are checks at points, not a proof.  Usage: python check_points.py umax npts [seed]"""
import sys, os, time, math, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from gzexact import Gz_exact, mvec, ldl_pivots

Q = flint.fmpq


def rat(x, den=1000):
    """a positive rational close to x (x > 0), with a few significant digits."""
    e = math.floor(math.log10(x))
    m = round(x / 10 ** e * den)
    return Q(m * 10 ** e, den) if e >= 0 else Q(m, den * 10 ** (-e))


def sample(u, it, rng):
    fam = it % 9
    wt = rat(math.exp(rng.uniform(-3, 3)))
    wU = [math.exp(rng.uniform(-3, 3)) for _ in range(u)]
    if fam == 7 and u >= 2:   # two U-weights tiny
        for h in rng.sample(range(u), 2):
            wU[h] *= 10 ** (-rng.uniform(2, 4))
    elif fam == 8:    # geometric spread of the U-weights
        r = math.exp(rng.uniform(0.2, 1.0))
        wU = [r ** h for h in range(u)]
    if fam == 1:      # one U-weight tiny
        wU[rng.randrange(u)] *= 10 ** (-rng.uniform(2, 4))
    elif fam == 2:    # all U-weights shrunk
        f = 10 ** (-rng.uniform(1, 4))
        wU = [x * f for x in wU]
    elif fam == 3:    # all U-weights large
        f = 10 ** (rng.uniform(1, 4))
        wU = [x * f for x in wU]
    elif fam == 4:    # one U-weight dominant
        wU[rng.randrange(u)] *= 10 ** (rng.uniform(1, 3))
    elif fam == 5:    # equal U-weights
        wU = [wU[0]] * u
    elif fam == 6:    # s_t tiny
        wt = rat(10 ** (-rng.uniform(1, 4)))
    return wt, [rat(x) for x in wU]


def check(wt, wU, gauge_check=False):
    u = len(wU)
    Gz = Gz_exact(wt, wU)
    st = wt * wt
    S = sum((x * x for x in wU), Q(0))
    v = mvec(st, S, u)
    n = u + 1
    Z = all(Gz[i, j] < 0 for i in range(n) for j in range(n) if i != j)
    Gv = [sum((Gz[i, j] * v[j] for j in range(n)), Q(0)) for i in range(n)]
    D = all(x > 0 for x in Gv)
    P = all(p > 0 for p in ldl_pivots(Gz))
    ratio = min(float(Gv[i] / (Gz[i, i] * v[i])) for i in range(n))
    g = True
    if gauge_check:
        g = (Gz_exact(wt, wU, gauge='r0') == Gz)
    return Z, D, P, g, ratio


if __name__ == '__main__':
    umax = int(sys.argv[1])
    npts = int(sys.argv[2])
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 2026
    rng = random.Random(seed)
    allok = True
    for u in range(1, umax + 1):
        t0 = time.time()
        nZ = nD = nP = nG = nGc = 0
        worst = 1.0
        for it in range(npts):
            wt, wU = sample(u, it, rng)
            Z, D, P, g, r = check(wt, wU, gauge_check=(it % 7 == 0))
            nZ += Z
            nD += D
            nP += P
            nG += g
            nGc += (it % 7 == 0)
            worst = min(worst, r)
            if not (Z and D and P and g):
                print(f'   FAILURE at u = {u}: w_t = {wt}, w_U = {wU}: Z {Z} D {D} P {P} gauge {g}', flush=True)
        ok = nZ == nD == nP == nG == npts
        allok &= ok
        print(f'u = {u}: {npts} rational points: (Z) {nZ}, (D) {nD}, pivots > 0 {nP}, gauge P_{{r0,h}} gives the same G_z at {nGc - (npts - nG)} of {nGc} checked points; '
              f'min_i (G_z v)_i/(G_ii v_i) = {worst:.3e}  ({time.time() - t0:.1f} s)  {"ok" if ok else "FAIL"}',
              flush=True)
    print('PASS' if allok else 'FAIL')
