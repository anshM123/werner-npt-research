"""Check, from our own objects, of the energy E_1 = A - B - 8 R_rank of NOTES.md and its link to the relaxed form Psi
(psi_exact.py, Lemma R1), in exact arithmetic, built from OUR objects:
 (1) R/2 = Pi diag(x V) Pi - sum_j (w P_j)(w Q_j)^T exactly, with V, P_j, Q_j transcribed from the
     measure-model formulas of relax_comparison.js (A_i = a_i/w_i, C_i = c_i/w_i, Z_i = zeta_i/w_i, free mass sum x_i = 1, mu_i = x_i);
 (2) the rational variance identity (5/18) B - N_D = 2 Var((x-b)V) + |E((x-2b)V)|^2 + d(E V^2 + (E V)^2), b = 13/72,
     d = 11/2592, and ||R/2||^2 <= (29/72)(B + 8 R_rank);
 (3) E_1 = A - B - 8 R_rank (A = our base form) is positive definite on {a perp w} at every sampled weight vector.
python bridge.py SEED NPER"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'typeI_G_m1_relax'))

import sys, random
import flint
from psi_exact import PsiExact, ldl_pivots
from test_A30 import parts

Q = flint.fmpq
Mt = flint.fmpq_mat


def sphere_point(rng, u):
    """rational w with sum w^2 = 1 and all w > 0 (inverse stereographic projection, retried until positive)"""
    while True:
        t = [Q(rng.randint(1, 30), rng.randint(31, 90)) for _ in range(u - 1)]
        n2 = sum((x * x for x in t), Q(0))
        w = [2 * x / (1 + n2) for x in t] + [(1 - n2) / (1 + n2)]
        if all(x > 0 for x in w):
            return w


def avg(mu, v):
    return sum((m * z for m, z in zip(mu, v)), Q(0))


def fields(ps, vec):
    """V, P_j, Q_j (the measure-model formulas, m = 1) at the amplitude vector vec = (alpha, c, rho, zeta) of psi_exact"""
    u, w = ps.u, ps.w
    na = u - 1
    alpha = vec[:na]; c = vec[na:na + u]; rho = vec[na + u]; zeta = vec[na + u + 1:]
    a = [sum((ps.Bm[i, r] * alpha[r] for r in range(na)), Q(0)) for i in range(u)]
    x = [wi * wi for wi in w]; mu = x
    A = [a[i] / w[i] for i in range(u)]; C = [c[i] / w[i] for i in range(u)]; Z = [zeta[i] / w[i] for i in range(u)]
    wt = ps.wt; sigma = wt * wt; K = sigma + Q(1, 2)
    p = avg(mu, C); F = avg(mu, [1 / (K + s) for s in x]); G = 2 + F
    V = [sigma * ((A[i] + C[i]) / (K + x[i]) - p / (2 * K * (K + x[i]))) for i in range(u)]
    st = sigma; k = K - st / 2
    b = st * rho + wt * avg(mu, Z)
    V = [V[i] + wt * Z[i] / (2 * (k + x[i])) - b / (4 * k * (k + x[i])) for i in range(u)]
    L = 2 + st / (k + st) + avg(mu, [1 / (k + s) for s in x])
    PP = [[sigma / (G * (K + s)) for s in x],
          [sigma / G * (x[i] * C[i] / (K + x[i]) + p / (2 * (K + x[i]))) for i in range(u)],
          [(wt * x[i] * Z[i] / (2 * (k + x[i])) + b / (4 * (k + x[i]))) / L for i in range(u)]]
    QQ = [[x[i] * A[i] / (K + x[i]) for i in range(u)], [1 / (K + s) for s in x], [1 / (k + s) for s in x]]
    cen = lambda v: [z - avg(mu, v) for z in v]
    P = [cen(v) for v in PP]; Qv = [cen(v) for v in QQ]
    return x, mu, V, P, Qv


def run(seed, nper):
    rng = random.Random(seed)
    ok = True
    for u in [2, 3, 4, 5, 6]:
        cnt = [0, 0, 0, 0]
        for it in range(nper):
            w = sphere_point(rng, u)
            wt = Q(rng.randint(1, 60), rng.randint(1, 60)) if it % 3 else Q(1, rng.randint(20, 200)) if it % 2 else Q(rng.randint(20, 200), 1)
            ps = PsiExact(wt, w)
            Amat, corr = parts(ps)
            n = Amat.nrows()
            na = u - 1
            # columns of R/2 for each basis vector of the amplitudes (linear in the amplitudes)
            def Rhalf(vec):
                alpha = vec[:na]; c = vec[na:na + u]; rho = vec[na + u]; zeta = vec[na + u + 1:]
                a = [sum((ps.Bm[i, r] * alpha[r] for r in range(na)), Q(0)) for i in range(u)]
                R = ps.Pi * (ps.Gcal(a, c) * ps.sig + ps.Hfull(rho, zeta)) * ps.Pi
                return [[R[r, h] / 2 for h in range(u)] for r in range(u)]
            # (1) decomposition identity at random amplitude vectors
            good1 = True
            for trial in range(3):
                vec = [Q(rng.randint(-9, 9), rng.randint(1, 7)) for _ in range(n)]
                x, mu, V, P, Qv = fields(ps, vec)
                Rh = Rhalf(vec)
                D = Mt(u, u, [x[i] * V[i] if i == j else Q(0) for i in range(u) for j in range(u)])
                PDP = ps.Pi * D * ps.Pi
                Kmat = Mt(u, u, [sum((w[i] * P[j][i] * w[l] * Qv[j][l] for j in range(3)), Q(0)) for i in range(u) for l in range(u)])
                diff = [PDP[i, l] - Kmat[i, l] - Rh[i][l] for i in range(u) for l in range(u)]
                good1 &= all(d == 0 for d in diff)
                # (2) variance identity and the source bound
                EV = avg(mu, V); B = avg(mu, [(x[i] + Q(1, 4)) * V[i] ** 2 for i in range(u)]) + EV * EV / 4
                ND = avg(mu, [x[i] * (1 - 2 * x[i]) * V[i] ** 2 for i in range(u)]) + avg(mu, [x[i] * V[i] for i in range(u)]) ** 2
                bb = Q(13, 72); dd = Q(11, 2592)
                y = [(x[i] - bb) * V[i] for i in range(u)]
                var = avg(mu, [z * z for z in y]) - avg(mu, y) ** 2
                rhs = 2 * var + avg(mu, [(x[i] - 2 * bb) * V[i] for i in range(u)]) ** 2 + dd * (avg(mu, [z * z for z in V]) + EV * EV)
                good1 &= (Q(5, 18) * B - ND == rhs)
                Rrank = sum((Kmat[i, l] ** 2 for i in range(u) for l in range(u)), Q(0))
                Rh2 = sum((Rh[i][l] ** 2 for i in range(u) for l in range(u)), Q(0))
                good1 &= Rh2 <= Q(29, 72) * (B + 8 * Rrank)
            cnt[0] += good1
            # (3) E_1 = A - B - 8 R_rank as a quadratic form in the amplitudes (polarisation on the basis)
            basis = [[Q(1) if i == j else Q(0) for i in range(n)] for j in range(n)]
            Vs = []; Ks = []
            for e in basis:
                x, mu, V, P, Qv = fields(ps, e)
                Vs.append(V)
                Ks.append([sum((w[i] * P[j][i] * w[l] * Qv[j][l] for j in range(3)), Q(0)) for i in range(u) for l in range(u)])
            Bm_ = [[avg(mu, [(x[i] + Q(1, 4)) * Vs[p][i] * Vs[q][i] for i in range(u)]) + avg(mu, Vs[p]) * avg(mu, Vs[q]) / 4
                    for q in range(n)] for p in range(n)]
            Rm_ = [[sum((Ks[p][r] * Ks[q][r] for r in range(u * u)), Q(0)) for q in range(n)] for p in range(n)]
            E1 = Mt(n, n, [Amat[p, q] - Bm_[p][q] - 8 * Rm_[p][q] for p in range(n) for q in range(n)])
            piv = ldl_pivots(E1)
            cnt[1] += all(z > 0 for z in piv)
            # and Psi - A/30 - (29/30) E_1 >= 0 (the chain), checked as a form
            T = Mt(n, n, [Amat[p, q] - corr[p, q] - Amat[p, q] / 30 - Q(29, 30) * E1[p, q] for p in range(n) for q in range(n)])
            # T = (29/30)(B + 8R) - corr should be PSD: test pivots >= 0 allowing zeros via a tiny check on eigen-sign
            pv = ldl_pivots(T + Mt(n, n, [Q(1, 10 ** 12) if p == q else Q(0) for p in range(n) for q in range(n)]))
            cnt[2] += all(z > 0 for z in pv)
            cnt[3] += 1
        print(f'u = {u}: decomposition + variance identity + source bound {cnt[0]}/{cnt[3]}; E_1 positive definite '
              f'{cnt[1]}/{cnt[3]}; (29/30)(B + 8R) - loss >= 0 {cnt[2]}/{cnt[3]}', flush=True)
        ok &= cnt[0] == cnt[1] == cnt[2] == cnt[3]
    print('PASS' if ok else 'FAIL')


if __name__ == '__main__':
    run(int(sys.argv[1]), int(sys.argv[2]))
