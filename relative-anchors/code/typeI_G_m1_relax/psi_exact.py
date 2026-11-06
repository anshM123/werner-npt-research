"""Exact (python-flint fmpq) construction of the relaxed form Psi of NOTES.md, for m = 1 (T = {t}), at rational
weights w_t, w_h (so all squared weights are rational squares).

  Psi(a, c, rho, zeta) = sig f(a, c) + <Hfull(rho, zeta), a w^T + w c^T> + Acal(rho, zeta)
                         - 1/4 sum_{r,h} R_rh^2 / lam_rh,     R = Pi (sig Gcal(a, c) + Hfull(rho, zeta)) Pi,
  a orthogonal to w (a = Bm alpha, Bm a rational basis of w-perp), Pi = I - w w^T/S,
  lam_rh = 5/12 (r != h), lam_hh = 3/2 - rho_h  (Lemma H of Theorem C).
Variables in this order: alpha (u - 1), c (u), rho, zeta (u).  The letter-2 coordinates are the last u + 1.
G_rel = Schur complement of Psi onto (rho, zeta); G_rel <= G (NOTES.md, Lemma R1)."""
import flint

Q = flint.fmpq
M = flint.fmpq_mat


def ident(n):
    return M(n, n, [Q(1) if i == j else Q(0) for i in range(n) for j in range(n)])


def diag(v):
    n = len(v)
    return M(n, n, [v[i] if i == j else Q(0) for i in range(n) for j in range(n)])


def col(v):
    return M(len(v), 1, list(v))


def outer(x, y):
    return M(len(x), len(y), [a * b for a in x for b in y])


class PsiExact:
    def __init__(self, wt, wU):
        self.wt = Q(wt)
        self.w = [Q(x) for x in wU]
        u = self.u = len(wU)
        w = self.w
        self.s = [x * x for x in w]
        sig = self.sig = self.wt * self.wt
        S = self.S = sum(self.s, Q(0))
        self.K = sig + S / 2
        self.kappa = (2 * sig + S) / (sig + S)
        self.d = [x + S / 2 + sig for x in self.s]
        E = diag(self.d) + outer(w, w) * Q(1, 2)
        self.Ei = E.inv()
        wtil = [self.wt] + w
        d0 = [sig + (sig + S) / 2] + [x + (sig + S) / 2 for x in self.s]
        E0 = diag(d0) + outer(wtil, wtil) * Q(1, 2)
        self.W = E0.inv()
        self.Pi = ident(u) - outer(w, w) * (1 / S)
        self.z = self.Ei * col(w)
        self.rho_ = [self.s[h] / self.d[h] for h in range(u)]
        self.lam = [[(Q(3, 2) - self.rho_[r]) if r == h else Q(5, 12) for h in range(u)] for r in range(u)]
        # basis of w-perp: columns w_0 e_r - w_r e_0, r = 1..u-1
        self.Bm = M(u, u - 1, [(w[0] if i == r else (-w[r] if i == 0 else Q(0))) for i in range(u) for r in range(1, u)])

    def Hfull(self, rho, zeta):
        u, w, wt, W = self.u, self.w, self.wt, self.W
        rho = Q(rho)
        wz = sum((w[h] * zeta[h] for h in range(u)), Q(0))
        Y0 = wt * rho + wz
        m = self.kappa * rho + wt * wz / (self.sig + self.S)
        c = [Y0 / 2 + w[r] * zeta[r] for r in range(u)]
        y = col([rho] + list(zeta))
        Wy = W * y
        H = ident(u) * (rho - m) + outer(w, w) * (m / (2 * self.K))
        H += M(u, u, [-c[r] * W[0, 1 + r] * w[h] + wt * c[r] * W[1 + h, 1 + r] for r in range(u) for h in range(u)])
        H += outer(w, w) * (-Wy[0, 0] / 2) + outer(w, [Wy[1 + h, 0] for h in range(u)]) * (wt / 2)
        return H

    def Gcal(self, a, c):
        u, w, Ei = self.u, self.w, self.Ei
        p = sum((w[h] * c[h] for h in range(u)), Q(0))
        return Ei * diag([2 * a[h] * w[h] for h in range(u)]) + diag([2 * c[h] * w[h] for h in range(u)]) * Ei \
            + (Ei - ident(u) * (1 / self.K)) * p

    def matrix(self):
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
        put(Faa, 0, 0, sig)
        put(Fac, 0, na, sig)
        put(Fac.transpose(), na, 0, sig)
        put(Mcc, na, na, sig)
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
                P[i][iy + j] += ga[i, 0] / 2
                P[iy + j][i] += ga[i, 0] / 2
            for i in range(u):
                P[na + i][iy + j] += gc[i, 0] / 2
                P[iy + j][na + i] += gc[i, 0] / 2
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
        Dl = diag([1 / self.lam[r][h] for r in range(u) for h in range(u)])
        corr = (L.transpose() * Dl * L) * Q(1, 4)
        Pm = M(n, n, [P[i][j] - corr[i, j] for i in range(n) for j in range(n)])
        return Pm, iy


def ldl_pivots(A):
    n = A.nrows()
    B = [[A[i, j] for j in range(n)] for i in range(n)]
    piv = []
    for p in range(n):
        pv = B[p][p]
        piv.append(pv)
        if pv == 0:
            return piv
        for i in range(p + 1, n):
            f = B[i][p] / pv
            if f != 0:
                Bi = B[i]
                Bp = B[p]
                for j in range(p + 1, n):
                    Bi[j] -= f * Bp[j]
    return piv
