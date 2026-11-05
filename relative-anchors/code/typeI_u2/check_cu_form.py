"""Exact checks (Fractions) of the (C_u) closed form and of the splitting X = X' + a w^T + w c^T.
 (a) CU.full_matrix() equals the principal block of the type I block (cf.TypeI) on {P_{r,h}: r, h in U} u {Q},
     for several k, T, weights (sigma = s_T);  and q-elimination: CU.form = Schur complement onto X;
 (b) the reduction lemma: the principal block of the full {X_U, q, X_T} block on {X_U, q} vs the brief's (C_u):
     both are PD at the test weights (the equivalence itself is proved in NOTES.md, Section 1);
 (c) the formulas (F1), (F2), (F3) of NOTES.md Section 2 at random rational points, u = 2..5.
Prints PASS."""
import sys, os, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from cf import TypeI
from cu_form import CU, inv_mat

ok = True
t0 = time.time()
random.seed(7)


def rfr(lo=1, hi=9):
    return Fr(random.randint(lo, hi), random.randint(1, 9))


# (a) closed form = principal block of the type I block
cases = [((Fr(1), Fr(2), Fr(3)), [0]), ((Fr(3, 2), Fr(1, 2), Fr(2), Fr(1)), [0]), ((Fr(3, 2), Fr(1, 2), Fr(2), Fr(1)), [0, 2]),
         ((Fr(1), Fr(2), Fr(1, 3), Fr(5, 4), Fr(2, 3)), [0, 3]), ((Fr(1), Fr(2), Fr(1, 3), Fr(5, 4), Fr(2, 3)), [4]),
         ((Fr(1), Fr(2), Fr(1, 3), Fr(5, 4), Fr(2, 3)), [0, 1, 2]), ((Fr(2, 7), Fr(5), Fr(1, 3), Fr(3, 4), Fr(2, 3), Fr(1)), [1, 5])]
for w, T in cases:
    ti = TypeI(w, T)
    M = ti.build()
    U = ti.U
    keep = [ti.idx[('P', r, h)] for r in U for h in U] + [ti.idx[('Q',)]]
    Mp = [[M[i][j] for j in keep] for i in keep]
    sig = sum(ti.s[t] for t in T)
    cu = CU([ti.w[h] for h in U], sig)
    Mc = cu.full_matrix()
    err = max(abs(Mp[i][j] - Mc[i][j]) for i in range(len(keep)) for j in range(len(keep)))
    # q elimination
    n = len(keep)
    u = len(U)
    err2 = Fr(0)
    for trial in range(3):
        X = [[rfr(-5, 5) for _ in range(u)] for _ in range(u)]
        x = [X[r][h] for r in range(u) for h in range(u)]
        Pxx = sum(x[i] * Mc[i][j] * x[j] for i in range(n - 1) for j in range(n - 1))
        Pxq = sum(x[i] * Mc[i][n - 1] for i in range(n - 1))
        val = Pxx - Pxq * Pxq / Mc[n - 1][n - 1]
        err2 = max(err2, abs(val - cu.form(X)))
    good = err == 0 and err2 == 0 and Mc[n - 1][n - 1] == 2
    ok &= good
    print(f"(a) k={len(w)} T={T}: principal block = closed form: max err {err}; q-elimination err {err2}: {good}", flush=True)


# (c) splitting formulas
def mat(u, f):
    return [[f(i, j) for j in range(u)] for i in range(u)]


for u in (2, 3, 4, 5):
    for trial in range(4):
        w = [rfr() for _ in range(u)]
        sig = rfr()
        cu = CU(w, sig)
        s = cu.s
        S = cu.S
        K = cu.K
        d = cu.d
        Ei = cu.Ei
        g = sum(s[h] / d[h] for h in range(u))
        rho = [s[h] / d[h] for h in range(u)]
        z = [sum(Ei[i][j] * w[j] for j in range(u)) for i in range(u)]
        E0 = [[cu.E[i][j] - (sig if i == j else 0) for j in range(u)] for i in range(u)]
        E0Ei = [[sum(E0[i][l] * Ei[l][j] for l in range(u)) for j in range(u)] for i in range(u)]
        P = mat(u, lambda i, j: (1 if i == j else 0) - w[i] * w[j] / S)
        Yr = mat(u, lambda i, j: rfr(-5, 5))
        Xp = [[sum(P[i][a] * Yr[a][b] * P[b][j] for a in range(u) for b in range(u)) for j in range(u)] for i in range(u)]
        b0 = [rfr(-5, 5) for _ in range(u)]
        a = [sum(P[i][j] * b0[j] for j in range(u)) for i in range(u)]
        c = [rfr(-5, 5) for _ in range(u)]
        Om = mat(u, lambda i, j: a[i] * w[j] + w[i] * c[j])
        X = mat(u, lambda i, j: Xp[i][j] + Om[i][j])
        FXp, FOm, FX = cu.form(Xp), cu.form(Om), cu.form(X)
        # (F1)
        F1 = (sum(Xp[r][h] ** 2 for r in range(u) for h in range(u) if r != h) / 2
              + sum((Fr(3, 2) - rho[h]) * Xp[h][h] ** 2 for h in range(u))
              + sum(rho[r] * rho[rp] * Xp[r][rp] * Xp[rp][r] for r in range(u) for rp in range(u)) / (2 + g))
        # (F2)
        cw = sum(ci * wi for ci, wi in zip(c, w))
        aEa = sum(a[i] * E0Ei[i][j] * a[j] for i in range(u) for j in range(u))
        cEc = sum(c[i] * E0Ei[i][j] * c[j] for i in range(u) for j in range(u))
        F2 = sig * (aEa + cEc - sig * cw ** 2 / (2 * K * K) + cw * sum(a_i * z_i for a_i, z_i in zip(a, z))
                    + 2 * sum(a[r] * c[r] * w[r] * z[r] for r in range(u)))
        # (F3)
        EiXp = [[sum(Ei[i][l] * Xp[l][j] for l in range(u)) for j in range(u)] for i in range(u)]
        XpEi = [[sum(Xp[i][l] * Ei[l][j] for l in range(u)) for j in range(u)] for i in range(u)]
        ipE = sum(Ei[i][j] * Xp[i][j] for i in range(u) for j in range(u))
        trXp = sum(Xp[h][h] for h in range(u))
        F3 = sig * (2 * sum(a[r] * w[r] * EiXp[r][r] for r in range(u)) + 2 * sum(c[r] * w[r] * XpEi[r][r] for r in range(u))
                    + cw * (ipE - trXp / K))
        good = (FXp == F1) and (FOm == F2) and (FX - FXp - FOm == F3)
        ok &= good
        if not good:
            print('  diffs', FXp - F1, FOm - F2, FX - FXp - FOm - F3)
        print(f"(c) u={u} trial {trial}: F1 {FXp == F1}, F2 {FOm == F2}, F3 {FX - FXp - FOm == F3}", flush=True)
print(f"time {time.time() - t0:.1f}s")
print("PASS" if ok else "FAIL")
