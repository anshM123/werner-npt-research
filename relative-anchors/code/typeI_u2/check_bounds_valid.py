"""Exact checks (Fractions) that the per-coordinate bounds of rowbounds2.py are valid at actual configurations, and
of the auxiliary inequalities used in NOTES.md Section 4.  Random rational normalised weights x_h (sum S), t = 1 - S/2,
u = 2..7, including spread and one-dominant configurations.  For every coordinate h:
  - the parameters (xi, theta, nu) lie in [0, 1];
  - pi_hk <= pair_Q(x_h, x_k, g) <= 0 for all k != h (needed for the signed Gershgorin row);
  - LB (bound) <= [Mc_hh + t sum_{k!=h} x_k pi_hk]/(1 - t) (exact signed weighted row of Mc);
  - max(ROWp, ROWm) >= weighted row sum of |Mca| and max(COLp, COLm) >= weighted column sum, both / (1 - t);
  - 0 <= (t/2) e_k m_k / Lam_k <= rho_k / G (used for |sigma_hk|) and the diagonal coupling of Ma <= (49/60) t rho_h.
Prints PASS."""
import sys, os, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from phi1 import blocks_norm, BETA
from rowbounds2 import bounds, pair_Q

random.seed(5)
ok = True
t0 = time.time()
ncfg = 0
for u in (2, 3, 4, 5, 7):
    for trial in range(40):
        S = Fr(random.randint(1, 199), 100)
        if trial % 5 == 0:
            S = Fr(1, random.randint(50, 2000))
        raw = [Fr(random.randint(1, 50)) ** random.choice([1, 2, 3]) for _ in range(u)]
        if trial % 3 == 1:
            raw[0] *= 1000
        tot = sum(raw)
        x = [r / tot * S for r in raw]
        t = 1 - S / 2
        # pi and sigma need only products ww[h][k] for the matrix; for row sums use x_k * pi_hk (rational)
        ww = [[None] * u for _ in range(u)]
        for h in range(u):
            for k in range(u):
                ww[h][k] = Fr(1)          # placeholder: we read pi_hk = Mc_hk / (t ww_hk) with ww = 1
        Ma, Mc, Mca = blocks_norm(x, t, ww)
        rho = [y / (1 + y) for y in x]; g = sum(rho); G = 2 + g
        eta2 = [y / (1 + y) ** 2 for y in x]; N = sum(eta2)
        Lam = [Fr(3, 2) - r for r in rho]
        e = [1 / (1 + x[h]) - eta2[h] / G for h in range(u)]
        m = [rho[h] + eta2[h] / G for h in range(u)]
        for h in range(u):
            # Ma diagonal coupling bound
            cpl = t * (BETA * rho[h] ** 2 * (N - eta2[h]) / G ** 2 + x[h] * e[h] ** 2 / Lam[h])
            ok &= cpl <= Fr(49, 60) * t * rho[h]
            for k in range(u):
                if k != h:
                    pi = Mc[h][k] / t
                    q = pair_Q(x[h], x[k], g, one=Fr(1))
                    ok &= (pi <= q) and (q <= 0)
                    ok &= 0 <= t / 2 * e[k] * m[k] / Lam[k] <= rho[k] / G
            sgn_row = (Mc[h][h] + t * sum(x[k] * Mc[h][k] / t for k in range(u) if k != h)) / (1 - t)
            row = sum(abs(Mca[h][k]) * x[k] for k in range(u) if k != h) + abs(Mca[h][h])
            col = sum(abs(Mca[k][h]) * x[k] for k in range(u) if k != h) + abs(Mca[h][h])
            Sp = S - x[h]; rSp = Sp / (1 + Sp); gp = g - rho[h]
            theta = (gp - rSp) / (Sp - rSp) if Sp - rSp != 0 else Fr(0)
            Dp = Sp - gp
            Nub_p = gp - Dp / (1 + Sp)
            nu = (N - eta2[h]) / Nub_p if Nub_p != 0 else Fr(0)
            ok &= 0 <= theta <= 1 and 0 <= nu <= 1 and 0 < x[h] / S <= 1
            b = bounds(S, x[h] / S, theta, nu, one=Fr(1))
            ok &= b['LB'] <= sgn_row
            ok &= max(b['ROWp'], b['ROWm']) >= row / (1 - t)
            ok &= max(b['COLp'], b['COLm']) >= col / (1 - t)
            ok &= b['Mca_hh'] == Mca[h][h]
        ncfg += 1
print(f"{ncfg} configurations, all bounds valid: {ok} ({time.time() - t0:.1f} s)")
print("PASS" if ok else "FAIL")
