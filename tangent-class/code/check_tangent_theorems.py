"""Checks: (1) charge formula 2^{n-c} q(X_{-c}) = sum_rho [...] with the vacuum row; (2) Theorem B (range and co-range
in T_e): q_n >= 2^{1-n} sum_{s!=t}||Y_st||^2; (3) Theorem C: min over a of q equals K^U(F) - 1/2||P_F alpha F||^2,
and EP on T_n <=> ||P_F alpha F||^2 <= 2 K^U(F); (4) Corollary D (rank-one blocks)."""
import numpy as np, itertools
from fock_tangent import *
rng = np.random.default_rng(11); worst = 0.0; ok = True
def charge_piece(X, n, k):
    st = strings(n); ex = np.array([sum(x != 0 for x in s) for s in st])
    M = np.zeros_like(X); mask = (ex[:, None] - ex[None, :]) == k; M[mask] = X[mask]; return M
for n in (2, 3, 4):
    st = strings(n); idx = {s: i for i, s in enumerate(st)}; ex = [sum(x != 0 for x in s) for s in st]
    for trial in range(3):
        U = [rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)) for _ in range(n)]; a = rng.normal(size=2) + 1j * rng.normal(size=2)
        u = tangent_range(U, a, n); v = rng.normal(size=(2, 3 ** n)) + 1j * rng.normal(size=(2, 3 ** n))
        X = np.outer(u[0], v[0].conj()) + np.outer(u[1], v[1].conj())
        # (1) charge decomposition and per-string formula
        tot = sum(q_direct(charge_piece(X, n, k), n) for k in range(-n, 2)); worst = max(worst, abs(tot - q_direct(X, n)) / abs(tot))
        for c in range(0, n):
            Xc = charge_piece(X, n, -c); qc = q_direct(Xc, n); val = 0.0
            for rho in [s for s in st if sum(x != 0 for x in s) == c]:
                H = [s for s in range(n) if rho[s] == 0]
                Y = {}
                for s_ in H:
                    for t_ in H:
                        B = np.zeros((2, 2), complex)
                        for z in (1, 2):
                            for zp in (1, 2):
                                row = [0] * n; row[s_] = z; col = list(rho); col[t_] = zp
                                B[z - 1, zp - 1] = X[idx[tuple(row)], idx[tuple(col)]]
                        Y[s_, t_] = B
                tau = {s_: np.trace(Y[s_, s_]) for s_ in H}; T = sum(tau.values()); x = X[0, idx[rho]]
                val += (4 / (c + 1)) * sum(np.sum(abs(Y[s_, t_]) ** 2) for s_ in H for t_ in H if s_ != t_) + 2 * sum(np.sum(abs(Y[s_, s_]) ** 2) for s_ in H) \
                       + abs(T - x) ** 2 - 2 * sum(abs(tau[s_]) ** 2 for s_ in H)
            err = abs(2 ** (n - c) * qc - val) / max(1, abs(val)); worst = max(worst, err)
        # (2) Theorem B: co-range restricted to H_0 + H_1
        P01 = np.diag([1.0 if e <= 1 else 0.0 for e in ex]); XB = X @ P01
        off = 0.0
        for s_ in range(n):
            for t_ in range(n):
                if s_ == t_: continue
                for z in (1, 2):
                    for zp in (1, 2):
                        r_ = [0] * n; r_[s_] = z; c_ = [0] * n; c_[t_] = zp; off += abs(XB[idx[tuple(r_)], idx[tuple(c_)]]) ** 2
        qB = q_direct(XB, n); bound = 2 ** (1 - n) * off
        if qB < bound - 1e-9 * max(1, abs(qB)): ok = False
        print(f"n={n} trial {trial}: Theorem B q={qB:.6f} >= bound {bound:.6f}")
print("charge-formula worst rel. error", worst)
# Theorem B adversarial: minimise (q - bound)/||X||^2 over random rank-two X with range in T_e and co-range in H_0+H_1
from scipy.optimize import minimize
for n in (3, 4):
    st = strings(n); idx = {s: i for i, s in enumerate(st)}; ex = np.array([sum(x != 0 for x in s) for s in st]); sel = np.where(ex <= 1)[0]
    def f(x):
        z = x[:len(x)//2] + 1j * x[len(x)//2:]
        U = [z[4*s:4*s+4].reshape(2, 2) for s in range(n)]; a = z[4*n:4*n+2]; vv = z[4*n+2:].reshape(2, len(sel))
        v = np.zeros((2, 3 ** n), complex); v[:, sel] = vv
        u = tangent_range(U, a, n); X = np.outer(u[0], v[0].conj()) + np.outer(u[1], v[1].conj())
        off = sum(abs(X[idx[tuple([z1 if w == s_ else 0 for w in range(n)])], idx[tuple([z2 if w == t_ else 0 for w in range(n)])]]) ** 2
                  for s_ in range(n) for t_ in range(n) if s_ != t_ for z1 in (1, 2) for z2 in (1, 2))
        return (q_direct(X, n) - 2 ** (1 - n) * off) / np.sum(abs(X) ** 2)
    dim = 2 * (4 * n + 2 + 2 * len(sel)); best = min(minimize(f, rng.normal(size=dim), method='L-BFGS-B', options={'maxiter': 300}).fun for _ in range(4))
    print(f"Theorem B adversarial n={n}: min (q - bound)/||X||^2 = {best:.3e}")
    if best < -1e-9: ok = False
# (3) Theorem C: minimisation over a
for n in (3, 4):
    for trial in range(3):
        U = [rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)) for _ in range(n)]
        O = Ops(U, np.zeros(2)); D = O.D; F = rng.normal(size=2 * D) + 1j * rng.normal(size=2 * D)
        K = np.vdot(F, (O.AbU + O.Vac - O.T + 0.5 * O.al.conj().T @ O.al) @ F).real
        Fr = np.stack([F[:D], F[D:]], 1); Qm, _ = np.linalg.qr(Fr); aF = O.al @ F; proj = np.linalg.norm(Qm.conj().T @ aF) ** 2
        pred = K - 0.5 * proj
        # direct: q(a) = <F|AbU+Vac-T|F> + 1/2||alF - sqrt2 sum conj(a_r) F_r||^2 minimised over a by least squares
        coef, *_ = np.linalg.lstsq(np.sqrt(2) * Fr, aF, rcond=None); a_opt = np.conj(coef)
        O2 = Ops(U, a_opt); direct = np.vdot(F, O2.Q @ F).real
        worst = max(worst, abs(pred - direct) / abs(direct))
        print(f"Theorem C n={n} trial {trial}: K - 1/2||P_F alpha F||^2 = {pred:.6f}, q at optimal a = {direct:.6f}")
# (4) rank-one blocks: lambda_min of the tangent form >= 0, n = 3, 4
for n in (3, 4):
    mn = 1e9
    for trial in range(10):
        U = [np.outer(rng.normal(size=2) + 1j * rng.normal(size=2), rng.normal(size=2) + 1j * rng.normal(size=2)) for _ in range(n)]
        Ug, ag = gauge(U, rng.normal(size=2) + 1j * rng.normal(size=2)); O = Ops(Ug, ag)
        M = O.Q * np.outer(np.concatenate([O.Dh, O.Dh]), np.concatenate([O.Dh, O.Dh])); mn = min(mn, np.linalg.eigvalsh(M)[0])
    print(f"rank-one blocks n={n}: min lambda_min = {mn:.3e}")
    if mn < -1e-12: ok = False
print("worst identity error", worst)
assert worst < 1e-10 and ok; print("PASS")
