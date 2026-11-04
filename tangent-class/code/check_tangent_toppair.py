"""Top pair {n-1, n}: (1) embedding identity  P(F_{n-1},F_n) = ||F||^2 + ||<a|F_{n-1}||^2 + Vac(F_{n-1}) - 1/2 ||D' F||^2,
D' = Lagrange map of Theorem 3 for n+1 copies, extra copy 0 with block w'_{0,r} = -sqrt2 a_r |1>  (so (D'F)_{0t} = -(alpha F_n + sqrt2 <a|F_{n-1})_t);
(2) the bound from Theorem 3 at n+1 copies:  P >= Vac(F_{n-1}) - ||<a|F_n||^2 + sum_r ||Gamma'_r F[r]||^2 (basis a = |a| e_1).
Both are checked against the chain form (Theorem A) restricted to the two sectors."""
import numpy as np, itertools
from fock_tangent import *
rng = np.random.default_rng(3); worst = 0; ok = True
for n in (3, 4):
    st = strings(n); idx = {s: i for i, s in enumerate(st)}
    for trial in range(4):
        U, a = rand_tangent(rng, n)
        # rotate reference so that a = |a| e_1 (unitary change of range basis keeps the gauge)
        th = np.angle(a[0]); W = np.array([[a[0], -np.conj(a[1])], [a[1], np.conj(a[0])]]) / np.linalg.norm(a)  # columns: a-hat, a-hat-perp
        U = [u @ np.conj(W) for u in U]; a = W.conj().T @ a   # new basis M = conj(W): U -> U M, a -> M^T a
        assert abs(a[1]) < 1e-12; a = np.array([abs(a[0]), 0]); p = a[0] ** 2
        O = Ops(U, a); D = O.D
        sel = np.concatenate([O.sector_idx(n - 1), O.sector_idx(n)])
        F = np.zeros(2 * D, complex); F[sel] = rng.normal(size=len(sel)) + 1j * rng.normal(size=len(sel))
        pair = np.vdot(F, O.Q @ F).real                              # tangent form on co-range in sectors n-1, n
        mask_top = np.concatenate([O.exc, O.exc]) == n; Fn = np.where(mask_top, F, 0); Fm = F - Fn
        Vac = np.vdot(Fm, O.Vac @ Fm).real
        aFm = O.abra @ Fm; aFn = O.abra @ Fn
        # D' F: pairs (s,t) of copies 1..n -> two-hole strings: (beta_s F)_{hole t part}... realised as T-L identity:
        # ||D_n F_{n-1}||^2 = 2 T(F_{n-1}) - ||alpha F_{n-1}||^2 on sector n-1 (Lagrange form), plus the copy-0 pairs:
        Dn = 2 * np.vdot(Fm, O.T @ Fm).real - np.linalg.norm(O.al @ Fm) ** 2
        D0 = np.linalg.norm(O.al @ Fn + np.sqrt(2) * O.abra @ Fm) ** 2
        Dp = Dn + D0
        emb = np.linalg.norm(F) ** 2 + np.linalg.norm(aFm) ** 2 + Vac - 0.5 * Dp
        e1 = abs(emb - pair) / abs(pair); worst = max(worst, e1)
        # Gamma'_r F[r] = Gamma_r G_r + xi_r P_r (copy 0 excited), xi = -sqrt2 a; Gamma_r G = r-th reference component of alpha^* G
        G = [Fm[r * D:(r + 1) * D] for r in range(2)]; Pv = [Fn[r * D:(r + 1) * D] for r in range(2)]
        gam = []
        for r in range(2):
            Gr = np.zeros(2 * D, complex); Gr[r * D:(r + 1) * D] = G[r]   # embed G_r as reference-r... alpha^* acts on strings:
            cre = O.al.conj().T @ G[r]                                     # alpha^* G_r in C^2 (x) H_n
            gam.append(np.linalg.norm(cre[r * D:(r + 1) * D] - np.sqrt(2) * a[r] * Pv[r]) ** 2)
        bound = Vac - np.linalg.norm(aFn) ** 2 + sum(gam)
        if pair < bound - 1e-9 * abs(pair): ok = False
        print(f"n={n} t={trial}: pair form {pair:.6f}, embedding {emb:.6f}, Thm-3(n+1) lower bound {bound:.6f}, p={p:.3f}")
print("embedding identity worst rel err", worst)
assert worst < 1e-10 and ok; print("PASS")
