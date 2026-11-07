"""Exact checks (python-flint, rational amplitudes; no floating point in any verdict) at rational weight vectors of
  (OD)  ||Ghat_tt'||_op <= H sqrt(theta_t theta_t'),      in the frame metric N_t = diag(1, I - w_U w_U^T/(sigma + S)),
        i.e.  H^2 theta_t theta_t' N_t - G_tt' N_t'^{-1} G_t't is positive semidefinite (exact LDL pivots);
  (S1)  ||Rhat_raw,tt'||_op <= c1 sqrt(theta_t theta_t'),  R_raw = form after eliminating psi, q and the common mean p
        (X_U along w_U w_U^T) only;
  (S2)  Shat_free,tt <= c2 theta_t I,  S_free = R_raw - G (positive semidefinite Gram);
and the diagonal bound of Theorem DG, Ghat_tt - (2/3)(1 - theta_t) I > 0, as a guard on the normalisation.
G is computed by direct elimination of the whole block (cf.TypeI, gauge P_{r0,h} deleted); R_raw from the formula of
lemmaq_exact.formula (itself checked against the direct G).  The largest ratios are also printed in floating point
(information only).
Usage: python exact_offdiag.py seed npoints H c1 c2      (H, c1, c2 as p/q)"""
import sys, random
import numpy as np
import flint
from lemmaq_exact import formula, direct_G, zeros, mat, sub

Q = flint.fmpq


def psd_exact(A):
    """exact test: A symmetric rational is positive semidefinite?  (LDL with pivoting on zero pivots)"""
    n = A.nrows()
    M = [[A[i, j] for j in range(n)] for i in range(n)]
    act = list(range(n))
    while act:
        # pick the largest diagonal among active to pivot
        p = max(act, key=lambda i: M[i][i])
        pv = M[p][p]
        if pv < 0:
            return False
        if pv == 0:
            # then the whole row must vanish on the active set
            if any(M[p][j] != 0 for j in act):
                return False
            act.remove(p)
            continue
        act.remove(p)
        for i in act:
            f = M[i][p] / pv
            if f != 0:
                for j in act:
                    M[i][j] -= f * M[p][j]
    return True


def Nmat(info):
    u, d, U, w, sig, S = info['u'], info['d'], info['U'], info['w'], info['sig'], info['S']
    N = zeros(d, d)
    N[0, 0] = 1
    for i, h in enumerate(U):
        for j, g in enumerate(U):
            N[i + 1, j + 1] = (1 if i == j else 0) - w[h] * w[g] / (sig + S)
    return N


def blk(M, t, tp, d):
    return mat([[M[t * d + i, tp * d + j] for j in range(d)] for i in range(d)])


def raw_split(wT, wU):
    G, Ptil, dl, info = formula(wT, wU)
    C, Ct, u, U, w, d, m = info['C'], info['Ct'], info['u'], info['U'], info['w'], info['d'], info['m']
    nz = C.nrows()
    V = zeros(nz, 2)
    V[nz - 1, 0] = 1
    for i, h in enumerate(U):
        for j, g in enumerate(U):
            V[i * u + j, 1] = w[h] * w[g]
    Cs = V.transpose() * C * V
    R = {}
    for t in range(m):
        for tp in range(m):
            full = Ct[t].transpose() * C.solve(Ct[tp])
            part = (V.transpose() * Ct[t]).transpose() * Cs.solve(V.transpose() * Ct[tp])
            R[(t, tp)] = G[(t, tp)] + full - part
    return G, R, info


def fl(A):
    return np.array([[float(A[i, j]) for j in range(A.ncols())] for i in range(A.nrows())])


def opratio(A, Nt, Ntp, kap2):
    """float: ||N_t^{-1/2} A N_t'^{-1/2}|| / sqrt(kap2)"""
    a, n1, n2 = fl(A), fl(Nt), fl(Ntp)
    e1, v1 = np.linalg.eigh(n1); e2, v2 = np.linalg.eigh(n2)
    i1 = v1 @ np.diag(e1 ** -0.5) @ v1.T; i2 = v2 @ np.diag(e2 ** -0.5) @ v2.T
    return np.linalg.norm(i1 @ a @ i2, 2) / np.sqrt(float(kap2))


def frac(sx):
    p, q = sx.split('/') if '/' in sx else (sx, '1')
    return Q(int(p), int(q))


if __name__ == '__main__':
    seed = int(sys.argv[1]); npts = int(sys.argv[2])
    H, c1, c2 = frac(sys.argv[3]), frac(sys.argv[4]), frac(sys.argv[5])
    rng = random.Random(seed)
    fails = {'OD': 0, 'S1': 0, 'S2': 0, 'DG': 0}
    worst = {'OD': 0, 'S1': 0, 'S2': 0}
    nb = 0
    for it in range(npts):
        m = rng.randint(2, 4); u = rng.randint(2, 6)
        mode = it % 4
        def amp(lo, hi):
            e = rng.uniform(lo, hi)
            x = 2.0 ** e
            return Q(int(round(x * 64)) or 1, 64)
        if mode == 0:
            wT = [amp(-3, 3) for _ in range(m)]; wU = [amp(-3, 3) for _ in range(u)]
        elif mode == 1:   # letter-2 mass dominant
            wT = [amp(2, 5) for _ in range(m)]; wU = [amp(-3, 0) for _ in range(u)]
        elif mode == 2:   # free mass dominant
            wT = [amp(-3, 0) for _ in range(m)]; wU = [amp(1, 4) for _ in range(u)]
        else:             # two comparable letter-2 copies, spread free weights
            base = amp(1, 3); wT = [base, base * Q(rng.randint(5, 9), 7)] + [amp(-2, 2) for _ in range(m - 2)]
            wU = [amp(-3, 3) for _ in range(u)]
        G, R, info = raw_split(wT, wU)
        Gd, _, _ = direct_G(wT, wU)
        d, s, sig = info['d'], info['s'], info['sig']
        N = Nmat(info)
        Ninv = N.inv()
        for t in range(m):
            assert blk(Gd, t, t, d) == G[(t, t)]
            th = s[t] / sig
            # Theorem DG guard: G_tt - (2/3)(1 - theta) N  > 0  (strict: all pivots > 0 -> use psd of (G - (2/3)(1-th)N) and det != 0)
            Dg = G[(t, t)] - Q(2, 3) * (1 - th) * N
            if not psd_exact(Dg) or Dg.det() == 0:
                fails['DG'] += 1
            S = R[(t, t)] - G[(t, t)]
            if not psd_exact(c2 * th * N - S):
                fails['S2'] += 1
            # float information: lambda_max(N^-1/2 S N^-1/2)/theta
            sf, nf = fl(S), fl(N)
            e, v = np.linalg.eigh(nf); ih = v @ np.diag(e ** -0.5) @ v.T
            worst['S2'] = max(worst['S2'], np.linalg.eigvalsh(ih @ sf @ ih)[-1] / float(th))
            for tp in range(m):
                if tp == t:
                    continue
                assert blk(Gd, t, tp, d) == G[(t, tp)]
                nb += 1
                kap = s[t] * s[tp] / sig ** 2
                A = H * H * kap * N - G[(t, tp)] * Ninv * G[(tp, t)]
                if not psd_exact(A):
                    fails['OD'] += 1
                B = c1 * c1 * kap * N - R[(t, tp)] * Ninv * R[(tp, t)]
                if not psd_exact(B):
                    fails['S1'] += 1
                worst['OD'] = max(worst['OD'], opratio(G[(t, tp)], N, N, kap))
                worst['S1'] = max(worst['S1'], opratio(R[(t, tp)], N, N, kap))
        print(f'point {it} (m={m}, u={u}, mode {mode}): cumulative failures {fails}', flush=True)
    print(f'{nb} off-diagonal blocks at {npts} rational points; H = {H}, c1 = {c1}, c2 = {c2}')
    print(f'exact failures: {fails}')
    print('float information, largest ratios: ||Ghat_tt\'||/sqrt(th th\') = %.4f, ||Rhat_raw||/sqrt(th th\') = %.4f, '
          'lambda_max(Shat_free,tt)/theta_t = %.4f' % (worst['OD'], worst['S1'], worst['S2']))
