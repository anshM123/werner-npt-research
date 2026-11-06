"""Check of the bound Psi >= A/30 + (29/30) E_1 > A/30 of NOTES.md, Section 4, i.e.
(29/30) A - corr is positive definite on {a perp w} x R^U x R x R^U, where A is the q-eliminated base at Omega (the
form P of psi_exact.py before the relaxation loss) and corr = (1/4) L^T diag(1/lam) L is the relaxed hard-part loss.
Also the coarser relaxation with every lam = 5/12.  Exact rational weights; exact LDL pivots (python-flint).
python test_A30.py SEED NPER [CNUM CDEN]   (with CNUM/CDEN = c the test is c A - corr > 0; default c = 29/30;
c = 1/2 is the negative control)"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'typeI_G_m1_relax'))

import sys, random
import flint
from psi_exact import PsiExact, ident, outer, col, diag, ldl_pivots

Q = flint.fmpq
M = flint.fmpq_mat


def parts(ps, coarse=False):
    """(A, corr) with Psi = A - corr; copied from PsiExact.matrix with the loss returned separately."""
    self = ps
    u, w, sig = self.u, self.w, self.sig
    na = u - 1
    n = na + u + 1 + u
    I = ident(u)
    A0 = I - self.Ei * sig
    zv = [self.z[h, 0] for h in range(u)]
    Mac = outer(zv, w) * Q(1, 2) + diag([w[h] * zv[h] for h in range(u)])
    Mcc = A0 - outer(w, w) * (sig / (2 * self.K ** 2))
    Bm = self.Bm
    Faa = Bm.transpose() * A0 * Bm
    Fac = Bm.transpose() * Mac
    P = [[Q(0)] * n for _ in range(n)]

    def put(blk, i0, j0, scale):
        for i in range(blk.nrows()):
            for j in range(blk.ncols()):
                P[i0 + i][j0 + j] += scale * blk[i, j]
    put(Faa, 0, 0, sig); put(Fac, 0, na, sig); put(Fac.transpose(), na, 0, sig); put(Mcc, na, na, sig)
    iy = na + u
    wcol = col(w)
    Hs = []
    for j in range(u + 1):
        rho = Q(1) if j == 0 else Q(0)
        zeta = [Q(1) if (j >= 1 and h == j - 1) else Q(0) for h in range(u)]
        H = self.Hfull(rho, zeta)
        Hs.append(H)
        ga = Bm.transpose() * (H * wcol)
        gc = H.transpose() * wcol
        for i in range(na):
            P[i][iy + j] += ga[i, 0] / 2; P[iy + j][i] += ga[i, 0] / 2
        for i in range(u):
            P[na + i][iy + j] += gc[i, 0] / 2; P[iy + j][na + i] += gc[i, 0] / 2
    mv = [self.kappa] + [self.wt * w[h] / (sig + self.S) for h in range(u)]
    for i in range(u + 1):
        for j in range(u + 1):
            base = Q(2) if (i == j == 0) else (Q(1) if i == j else Q(0))
            P[iy + i][iy + j] += base - mv[i] * mv[j] / 2
    Pi = self.Pi
    Rcols = []
    for k in range(n):
        if k < na:
            a = [Bm[i, k] for i in range(u)]
            R = Pi * (self.Gcal(a, [Q(0)] * u) * sig) * Pi
        elif k < na + u:
            c = [Q(1) if i == k - na else Q(0) for i in range(u)]
            R = Pi * (self.Gcal([Q(0)] * u, c) * sig) * Pi
        else:
            R = Pi * Hs[k - iy] * Pi
        Rcols.append([R[r, h] for r in range(u) for h in range(u)])
    L = M(u * u, n, [Rcols[k][i] for i in range(u * u) for k in range(n)])
    lam = [[Q(5, 12) for h in range(u)] for r in range(u)] if coarse else self.lam
    Dl = diag([1 / lam[r][h] for r in range(u) for h in range(u)])
    corr = (L.transpose() * Dl * L) * Q(1, 4)
    A = M(n, n, [P[i][j] for i in range(n) for j in range(n)])
    return A, corr


def rq(rng, lo, hi):
    """a rational square root weight: w = p/q with w^2 spread over e^[lo, hi] roughly"""
    import math
    x = math.exp(rng.uniform(lo, hi))
    q = rng.randint(1, 40)
    p = max(1, round(x * q))
    return Q(p, q)


if __name__ == '__main__':
    seed, nper = int(sys.argv[1]), int(sys.argv[2])
    CNUM, CDEN = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (29, 30)
    rng = random.Random(seed)
    ok = True
    for u in [2, 3, 4, 5, 6, 8]:
        good = 0; tot = 0; worst = None
        for it in range(nper):
            mode = it % 4
            if mode == 0:
                wU = [rq(rng, -2, 2) for _ in range(u)]; wt = rq(rng, -2, 2)
            elif mode == 1:
                wU = [rq(rng, -4, 4) for _ in range(u)]; wt = rq(rng, -6, -3)      # small letter weight
            elif mode == 2:
                wU = [rq(rng, -4, 4) for _ in range(u)]; wt = rq(rng, 3, 6)        # large letter weight
            else:
                wU = [rq(rng, 0, 0.3) for _ in range(u)]; wU[0] = rq(rng, 3, 5); wt = rq(rng, -1, 1)   # one dominant free weight
            ps = PsiExact(wt, wU)
            A, corr = parts(ps)
            T1 = A * Q(CNUM, CDEN) - corr
            p1 = ldl_pivots(T1)
            A2, corr2 = parts(ps, coarse=True)
            T2 = A2 * Q(CNUM, CDEN) - corr2
            p2 = ldl_pivots(T2)
            r1 = all(x > 0 for x in p1); r2 = all(x > 0 for x in p2)
            tot += 1; good += r1 and r2
            if not (r1 and r2):
                print(f'  u={u} wt={wt} wU={[str(x) for x in wU]}: (29/30)A - corr PD {r1}; coarse PD {r2}', flush=True)
        ok &= good == tot
        print(f'u = {u}: ({CNUM}/{CDEN}) A - corr > 0 and the coarse version > 0 at {good}/{tot} weight vectors', flush=True)
    print('PASS' if ok else 'FAIL')
