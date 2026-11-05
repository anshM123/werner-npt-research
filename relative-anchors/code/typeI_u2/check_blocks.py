"""Exact checks (Fractions):
 (1) on a perp w, Phi1 (unnormalised: soft form (F2) minus the coupling bound built from (F3)) equals the quadratic
     form of the normalised block formulas [[Ma, Mca^T], [Mca, Mc]] (NOTES.md Section 3), u = 2..5, random rational
     weights and sigma, random rational (a perp w, c);
 (2) the chain at random points: F(X) - sigma Phi1(a, c) >= 0 for X = X' + a w^T + w c^T (exact, Fractions), where F is
     the (C_u) form of cu_form.py; and the coupling bound K^2 ||Gm||_F^2 <= 29 S (|a|^2 + |c|^2) of Corollary Q.
Prints PASS."""
import sys, os, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from cu_form import CU
from phi1 import phi1_unnorm, blocks_norm, coupling_matrix

random.seed(11)
ok = True
t0 = time.time()


def rfr(lo, hi):
    return Fr(random.randint(lo, hi), random.randint(1, 9))


for u in (2, 3, 4, 5):
    for trial in range(3):
        w = [rfr(1, 12) for _ in range(u)]
        sig = rfr(1, 40)
        s = [x * x for x in w]; S = sum(s); K = sig + S / 2
        x = [sh / K for sh in s]; t = sig / K
        ww = [[w[h] * w[k] / K for k in range(u)] for h in range(u)]
        Ma, Mc, Mca = blocks_norm(x, t, ww)
        P = [[(1 if i == j else 0) - w[i] * w[j] / S for j in range(u)] for i in range(u)]
        for rep in range(2):
            b0 = [rfr(-6, 6) for _ in range(u)]
            a = [sum(P[i][j] * b0[j] for j in range(u)) for i in range(u)]
            c = [rfr(-6, 6) for _ in range(u)]
            lhs = phi1_unnorm(w, sig, a, c)
            rhs = (sum(a[i] * Ma[i][j] * a[j] for i in range(u) for j in range(u))
                   + sum(c[i] * Mc[i][j] * c[j] for i in range(u) for j in range(u))
                   + 2 * sum(c[i] * Mca[i][j] * a[j] for i in range(u) for j in range(u)))
            good = lhs == rhs
            # (2) the chain
            cu = CU(w, sig)
            Yr = [[rfr(-5, 5) for _ in range(u)] for _ in range(u)]
            Xp = [[sum(P[i][p] * Yr[p][q] * P[q][j] for p in range(u) for q in range(u)) for j in range(u)] for i in range(u)]
            X = [[Xp[i][j] + a[i] * w[j] + w[i] * c[j] for j in range(u)] for i in range(u)]
            gap = cu.form(X) - sig * lhs
            good2 = gap >= 0
            # Corollary Q: K^2 ||Gm||_F^2 <= 29 S (|a|^2 + |c|^2)  (normalised: ||Gm|| <= sqrt(29) ||Omega||)
            Gm = coupling_matrix(w, sig, a, c)
            good3 = K * K * sum(v * v for row in Gm for v in row) <= 29 * S * (sum(v * v for v in a) + sum(v * v for v in c))
            good2 = good2 and good3
            ok &= good and good2
            print(f"u={u} trial {trial}.{rep}: Phi1 = block form: {good};  F(X) - sigma Phi1 = {float(gap):.4g} >= 0: {good2}", flush=True)
print(f"time {time.time() - t0:.1f}s")
print("PASS" if ok else "FAIL")
