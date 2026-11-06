"""Exploration (Psi built exactly by psi_exact.py, analysed in floating point; no claim rests on this).
Coordinates: global (gamma, xi, rho) with c = gamma w/|w| + c_perp, zeta = xi w/|w| + zeta_perp; perp (a, c_perp,
zeta_perp), each in an orthonormal basis of w-perp.  After Jacobi scaling, report lambda_min of Psi, of the global
3 x 3 block, of the perp block, the norm of the coupling, and the Schur ratio |C_gp|^2/(lam_glob lam_perp)."""
import sys, math, random
import numpy as np
import flint
from psi_exact import PsiExact
from certify_points import family

rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 3)
us = [int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else '4,6,8').split(',')]
for u in us:
    worst = dict(lmin=9, glob=9, perp=9, ratio=0)
    for k in range(16):
        wt, wU = family(u, k, rng)
        P = PsiExact(wt, wU)
        Mx, iy = P.matrix()
        n = Mx.nrows()
        M = np.array([[float(Mx[i, j]) for j in range(n)] for i in range(n)])
        w = np.array([float(x) for x in wU])
        Bm = np.array([[float(P.Bm[i, j]) for j in range(u - 1)] for i in range(u)])
        wh = w / np.linalg.norm(w)
        Q = np.linalg.svd(wh[None, :])[2][1:].T          # orthonormal basis of w-perp (u x (u-1))
        # psi_exact coordinates: alpha (a = Bm alpha), c (u), rho, zeta (u)
        # new coordinates: gamma, xi, rho, a_perp (u-1), c_perp (u-1), zeta_perp (u-1)
        T = np.zeros((n, 3 + 3 * (u - 1)))
        Bpinv = np.linalg.pinv(Bm)                        # alpha = Bpinv a for a in w-perp
        T[:u - 1, 3:3 + (u - 1)] = Bpinv @ Q
        T[u - 1:2 * u - 1, 0] = wh
        T[u - 1:2 * u - 1, 3 + (u - 1):3 + 2 * (u - 1)] = Q
        T[2 * u - 1, 2] = 1.0
        T[2 * u:, 1] = wh
        T[2 * u:, 3 + 2 * (u - 1):] = Q
        Mt = T.T @ M @ T
        dsc = 1 / np.sqrt(np.diag(Mt))
        Mn = Mt * np.outer(dsc, dsc)
        lmin = np.linalg.eigvalsh(Mn).min()
        Gb = Mn[:3, :3]; Pp = Mn[3:, 3:]; Cg = Mn[:3, 3:]
        lg = np.linalg.eigvalsh(Gb).min(); lp = np.linalg.eigvalsh(Pp).min()
        nc = np.linalg.norm(Cg, 2)
        ratio = nc ** 2 / (lg * lp) if lg > 0 and lp > 0 else float('inf')
        worst['lmin'] = min(worst['lmin'], lmin); worst['glob'] = min(worst['glob'], lg)
        worst['perp'] = min(worst['perp'], lp); worst['ratio'] = max(worst['ratio'], ratio)
    print(f"u = {u}: Jacobi-scaled lambda_min(Psi) >= {worst['lmin']:.3e}; global block >= {worst['glob']:.3e}; "
          f"perp block >= {worst['perp']:.3e}; max |C_gp|^2/(lam_g lam_p) = {worst['ratio']:.3e}", flush=True)
