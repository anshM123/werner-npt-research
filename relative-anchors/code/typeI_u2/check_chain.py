"""End-to-end exact test of steps (ii)-(iii) of Theorem C (NOTES.md Section 4): at random rational weights and sigma
(u = 2..6), with the actual per-coordinate quantities LB_h, ROW_h, COL_h (exact, no moment bounds) and with the bounded
ones of rowbounds2.py, the matrix
      Phi1 - (1 - t) diag( 1 - COL_k/tau (a-block),  LB_h - tau ROW_h (c-block) ),   tau = 4/5,
is positive semidefinite (exact LDL^T in Fractions), and the diagonal is positive.  Prints PASS."""
import sys, os, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from phi1 import blocks_norm
from rowbounds2 import bounds

random.seed(9)
TAU = Fr(4, 5)


def psd(M):
    """exact: symmetric M is PSD (LDL^T with zero-pivot handling)."""
    n = len(M)
    A = [row[:] for row in M]
    for i in range(n):
        if A[i][i] < 0:
            return False
        if A[i][i] == 0:
            if any(A[i][j] != 0 for j in range(i + 1, n)):
                return False
            continue
        for j in range(i + 1, n):
            f = A[j][i] / A[i][i]
            for k in range(i, n):
                A[j][k] -= f * A[i][k]
    return True


ok = True
t0 = time.time()
cnt = 0
for u in (2, 3, 4, 6):
    for trial in range(8):
        w = [Fr(random.randint(1, 30), random.randint(1, 9)) for _ in range(u)]
        if trial % 3 == 1:
            w[0] *= 20
        sig = Fr(random.randint(1, 400), random.randint(1, 40))
        s = [y * y for y in w]; S0 = sum(s); K = sig + S0 / 2
        x = [y / K for y in s]; t = sig / K; Sx = sum(x)
        ww = [[w[h] * w[k] / K for k in range(u)] for h in range(u)]
        Ma, Mc, Mca = blocks_norm(x, t, ww)
        n = 2 * u
        Phi = [[Fr(0)] * n for _ in range(n)]
        for i in range(u):
            for j in range(u):
                Phi[i][j] = Ma[i][j]; Phi[u + i][u + j] = Mc[i][j]
                Phi[u + i][j] = Mca[i][j]; Phi[j][u + i] = Mca[i][j]
        rho = [y / (1 + y) for y in x]; g = sum(rho); eta2 = [y / (1 + y) ** 2 for y in x]; N = sum(eta2)
        for mode in ('exact', 'bounded'):
            D = [Fr(0)] * n
            for h in range(u):
                # actual quantities: x_k pi_hk = Mc_hk * (w_k/w_h) / ... -> use ww: Mc_hk w_k/w_h = Mc_hk * x_k / ww_hk
                LB = (Mc[h][h] + sum(Mc[h][k] * x[k] / ww[h][k] for k in range(u) if k != h)) / (1 - t)
                ROW = (abs(Mca[h][h]) + sum(abs(Mca[h][k]) * x[k] / ww[h][k] for k in range(u) if k != h)) / (1 - t)
                COL = (abs(Mca[h][h]) + sum(abs(Mca[k][h]) * x[k] / ww[h][k] for k in range(u) if k != h)) / (1 - t)
                if mode == 'bounded':
                    Sp = Sx - x[h]; rSp = Sp / (1 + Sp); gp = g - rho[h]
                    th = (gp - rSp) / (Sp - rSp) if Sp != rSp else Fr(0)
                    Dp = Sp - gp; Nub_p = gp - Dp / (1 + Sp)
                    nu = (N - eta2[h]) / Nub_p if Nub_p != 0 else Fr(0)
                    bb = bounds(Sx, x[h] / Sx, th, nu, one=Fr(1))
                    ok &= bb['LB'] <= LB and max(bb['ROWp'], bb['ROWm']) >= ROW and max(bb['COLp'], bb['COLm']) >= COL
                    LB, ROW, COL = bb['LB'], max(bb['ROWp'], bb['ROWm']), max(bb['COLp'], bb['COLm'])
                D[h] = (1 - t) * (1 - COL / TAU)
                D[u + h] = (1 - t) * (LB - TAU * ROW)
            ok &= all(dd > 0 for dd in D)
            R = [[Phi[i][j] - (D[i] if i == j else 0) for j in range(n)] for i in range(n)]
            good = psd(R)
            ok &= good
            cnt += 1
        print(f"u={u} trial {trial}: Phi1 - (1-t) diag(...) PSD (exact and bounded quantities): {good}", flush=True)
print(f"{cnt} matrices, time {time.time() - t0:.1f}s")
print("PASS" if ok else "FAIL")
