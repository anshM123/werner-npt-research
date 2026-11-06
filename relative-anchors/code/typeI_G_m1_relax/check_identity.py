"""Exact check (python-flint fmpq) of the decomposition behind Lemma R1 (NOTES.md, Section 2), against the closed
form cf.TypeI of the whole type I block (m = 1, T = {0}, U = {1..u}, gauge x_t = 0, units of 2Y):

  (I1) min_q Phi(X, q, rho, zeta) = Phi_q(X, 0, 0) + <Hfull(rho, zeta), X> + Acal(rho, zeta)   for all X, rho, zeta;
  (I2) Phi_q(Omega, 0, 0) = sig f(a, c)   for Omega = a w^T + w c^T, a orthogonal to w   (Theorem C, (F2));
  (I3) Phi_q(X' + Omega, 0, 0) - Phi_q(X', 0, 0) - Phi_q(Omega, 0, 0) = sig <Gcal(a, c), X'>  for hard X'   ((F3));
  (I4) Phi_q(X', 0, 0) - sum lam_rh X'_rh^2 >= 0 on the hard space (Lemma H), tested by exact LDL^T pivots;
  (I5) Psi(a, c, rho, zeta) <= Phi_q(X' + Omega, rho, zeta) at random rational points (the inequality of Lemma R1);
  (I6) G - G_rel is positive semidefinite (exact LDL^T, symmetric pivoting not needed generically).
Random rational weights (w rational), random rational test vectors; u = 2..6."""
import sys, random
import flint
from cf import TypeI
from psi_exact import PsiExact, ldl_pivots, ident, outer, col, diag

Q = flint.fmpq
M = flint.fmpq_mat
random.seed(12345)


def rq(lo=1, hi=30):
    return Q(random.randint(lo, hi), random.randint(lo, hi))


def rv(n):
    return [Q(random.randint(-9, 9), random.randint(1, 7)) for _ in range(n)]


def full_block(wt, wU):
    w = [Q(wt)] + [Q(x) for x in wU]
    ti = TypeI(w, [0], F=Q)
    Mb = ti.build()
    U = ti.U
    el = [ti.idx[('P', r, h)] for r in U for h in U]
    iq = ti.idx[('Q',)]
    tb = [ti.idx[('R', 0)]] + [ti.idx[('Z', 0, h)] for h in U]
    return Mb, el, iq, tb


def phi_q(Mb, el, iq, tb, X, rho, zeta):
    """min over q of the block form at (X, q, rho, zeta) (exact; Phi_qq = 2 > 0)."""
    u = len(zeta)
    vec = {}
    for r in range(u):
        for h in range(u):
            vec[el[r * u + h]] = X[r][h]
    vec[tb[0]] = Q(rho)
    for h in range(u):
        vec[tb[1 + h]] = zeta[h]
    keys = list(vec)
    val = sum((vec[i] * Mb[i][j] * vec[j] for i in keys for j in keys), Q(0))
    lin = sum((Mb[iq][j] * vec[j] for j in keys), Q(0))
    return val - lin * lin / Mb[iq][iq]


if __name__ == '__main__':
    ok = True
    for u in range(2, 7):
        for trial in range(3):
            wt = rq(); wU = [rq() for _ in range(u)]
            Mb, el, iq, tb = full_block(wt, wU)
            P = PsiExact(wt, wU)
            w = P.w
            S = P.S
            Pi = P.Pi
            sig = P.sig
            # random hard X' = Pi Y Pi, a = Pi b, c, rho, zeta
            Y = M(u, u, [x for x in rv(u * u)])
            Xh = Pi * Y * Pi
            av = Pi * col(rv(u))
            a = [av[i, 0] for i in range(u)]
            assert sum((w[h] * a[h] for h in range(u)), Q(0)) == 0
            c = rv(u); rho = rv(1)[0]; zeta = rv(u)
            Om = outer(a, w) + outer(w, c)
            L = lambda A: [[A[i, j] for j in range(u)] for i in range(u)]
            z0 = [Q(0)] * u
            # (I1) on a random full X
            Xf = M(u, u, rv(u * u))
            H = P.Hfull(rho, zeta)
            lhs = phi_q(Mb, el, iq, tb, L(Xf), rho, zeta)
            Hx = sum((H[i, j] * Xf[i, j] for i in range(u) for j in range(u)), Q(0))
            rhs = phi_q(Mb, el, iq, tb, L(Xf), 0, z0) + Hx + (phi_q(Mb, el, iq, tb, L(Xf * 0), rho, zeta))
            mv = [P.kappa] + [P.wt * w[h] / (sig + S) for h in range(u)]
            mm = mv[0] * rho + sum((mv[1 + h] * zeta[h] for h in range(u)), Q(0))
            Acal = 2 * rho * rho + sum((x * x for x in zeta), Q(0)) - mm * mm / 2
            i1 = (lhs == rhs) and (phi_q(Mb, el, iq, tb, L(Xf * 0), rho, zeta) == Acal)
            # (I2)
            Ei = P.Ei
            A0 = ident(u) - Ei * sig
            zv = [P.z[h, 0] for h in range(u)]
            p = sum((w[h] * c[h] for h in range(u)), Q(0))
            za = sum((zv[h] * a[h] for h in range(u)), Q(0))
            f = (M(1, u, a) * A0 * col(a))[0, 0] + (M(1, u, c) * A0 * col(c))[0, 0] - sig * p * p / (2 * P.K ** 2) \
                + p * za + 2 * sum((a[h] * c[h] * w[h] * zv[h] for h in range(u)), Q(0))
            i2 = phi_q(Mb, el, iq, tb, L(Om), 0, z0) == sig * f
            # (I3)
            Gc = P.Gcal(a, c)
            lhs3 = phi_q(Mb, el, iq, tb, L(Xh + Om), 0, z0) - phi_q(Mb, el, iq, tb, L(Xh), 0, z0) - \
                phi_q(Mb, el, iq, tb, L(Om), 0, z0)
            i3 = lhs3 == sig * sum((Gc[i, j] * Xh[i, j] for i in range(u) for j in range(u)), Q(0))
            # (I4) F - lam-diagonal on the hard space: hard basis Bh (vec of Pi E_rh Pi for r,h < u-1 suffices to span)
            Bm = P.Bm
            hb = []
            for i in range(u - 1):
                for j in range(u - 1):
                    Eij = M(u, 1, [Bm[r, i] for r in range(u)]) * M(1, u, [Bm[r, j] for r in range(u)])
                    hb.append(Eij)
            nh = len(hb)
            Fm = M(nh, nh)
            def Fq(A):
                return phi_q(Mb, el, iq, tb, L(A), 0, z0) - sum((P.lam[r][h] * A[r, h] ** 2 for r in range(u)
                                                                  for h in range(u)), Q(0))
            dgs = [Fq(hb[i]) for i in range(nh)]
            for i in range(nh):
                Fm[i, i] = dgs[i]
                for j in range(i + 1, nh):
                    v = (Fq(hb[i] + hb[j]) - dgs[i] - dgs[j]) / 2
                    Fm[i, j] = v
                    Fm[j, i] = v
            pv = ldl_pivots(Fm)
            i4 = len(pv) == nh and all(x > 0 for x in pv)
            # (I5): Psi(a, c, rho, zeta) <= Phi_q(X' + Omega, rho, zeta)
            Pm, iy = P.matrix()
            # coordinates: alpha with a = Bm alpha.  Solve Bm alpha = a (a in w-perp): least squares exact via normal eqs
            BtB = Bm.transpose() * Bm
            alpha = BtB.solve(Bm.transpose() * col(a))
            vv = [alpha[i, 0] for i in range(u - 1)] + list(c) + [rho] + list(zeta)
            n = len(vv)
            psi_val = sum((vv[i] * Pm[i, j] * vv[j] for i in range(n) for j in range(n)), Q(0))
            phi_val = phi_q(Mb, el, iq, tb, L(Xh + Om), rho, zeta)
            i5 = psi_val <= phi_val
            # (I6) G - G_rel >= 0
            # exact G from the full block (gauge x_t = 0) and G_rel = Schur complement of Psi
            elq = el + [iq]
            Cm = M(len(elq), len(elq), [Mb[i][j] for i in elq for j in elq])
            Bc = M(len(elq), len(tb), [Mb[i][j] for i in elq for j in tb])
            Am = M(len(tb), len(tb), [Mb[i][j] for i in tb for j in tb])
            G = Am - Bc.transpose() * Cm.solve(Bc)
            na = iy
            A11 = M(na, na, [Pm[i, j] for i in range(na) for j in range(na)])
            A12 = M(na, u + 1, [Pm[i, iy + j] for i in range(na) for j in range(u + 1)])
            A22 = M(u + 1, u + 1, [Pm[iy + i, iy + j] for i in range(u + 1) for j in range(u + 1)])
            Grel = A22 - A12.transpose() * A11.solve(A12)
            pv6 = ldl_pivots(G - Grel)
            Dm = G - Grel
            if len(pv6) == u + 1 and all(x > 0 for x in pv6):
                i6 = True
            else:
                # positive semidefinite test with a possible kernel: all principal minors of order <= u + 1 >= 0
                from itertools import combinations
                i6 = True
                for kk in range(1, u + 2):
                    for idx in combinations(range(u + 1), kk):
                        sub = M(kk, kk, [Dm[i, j] for i in idx for j in idx])
                        if sub.det() < 0:
                            i6 = False
                print('    G - G_rel singular: pivots', [float(x) for x in pv6], 'all principal minors >= 0:', i6)
            pvG = ldl_pivots(Grel)
            allok = i1 and i2 and i3 and i4 and i5 and i6
            ok &= allok
            print(f'u = {u}, trial {trial}: I1 {i1} I2 {i2} I3 {i3} I4 {i4} I5 {i5} I6 {i6}; G_rel pivots > 0: '
                  f'{len(pvG) == u + 1 and all(x > 0 for x in pvG)}', flush=True)
    print('PASS' if ok else 'FAIL')
