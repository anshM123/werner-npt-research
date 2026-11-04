"""Numerical checks for Proposition 10 (README): the two-copy anchor phi_u = |10> + u|01>, and the sector statement.

One copy: M[(a,y),(b,y')] = d_ab d_yy' - 1/2 d_ay d_by' (row a, column y), so q_n(X) = <vec X, M^{(x)n} vec X>.
For a row vector p: V_p = p (x) I, D_p = V_p^* M V_p, R_p = M V_p D_p^{-1} V_p^* M, A = M^Gamma (Gamma = transpose on
the row factor) and T_p = A^{-1/2} R_p^Gamma A^{-1/2}.
(a) T for the vacuum on one copy has spectrum {-1/3, 1}.
(b) T_{phi_u} lies in [-3, 1] on a grid of u (the exact proof is check_pair_anchor_exact.py), and D > 0.
(c) The same for random complex two-copy vectors x (x) |0> + |0> (x) y (local unitary invariance).
(d) Sector level: if A_0 + A_1 has the eigenvalue 1, every sector form is >= 0 (n = 5, 6, random data)."""
import numpy as np, itertools
from threelevel_certificate import GeneralModel, gauge

def Mel(a, y, b, yp):
    v = 1.0
    for t in range(len(a)):
        v *= (1.0 if (a[t] == b[t] and y[t] == yp[t]) else 0.0) - (0.5 if (a[t] == y[t] and b[t] == yp[t]) else 0.0)
        if v == 0: return 0.0
    return v

def Tspec(p, n):
    S = list(itertools.product(range(3), repeat=n)); N = len(S)
    idx = [(a, y) for a in S for y in S]; ip = {k: i for i, k in enumerate(idx)}
    M = np.array([[Mel(a, y, b, yp) for (b, yp) in idx] for (a, y) in idx])
    Vp = np.zeros((N * N, N), complex)
    for a, pa in p.items():
        for j, x in enumerate(S): Vp[ip[(a, x)], j] += pa
    D = Vp.conj().T @ M @ Vp
    R = M @ Vp @ np.linalg.inv(D) @ Vp.conj().T @ M
    perm = np.array([[ip[(b, yp)] for (b, y) in idx] for (a, yp) in idx])
    perm2 = np.array([[ip[(a, y)] for (b, y) in idx] for (a, yp) in idx])
    A = M[perm, perm2]; RG = R[perm, perm2]
    w, V = np.linalg.eigh(A); Ais = V @ np.diag(w ** -0.5) @ V.T
    T = Ais @ RG @ Ais; T = (T + T.conj().T) / 2
    ev = np.linalg.eigvalsh(T)
    return ev[0], ev[-1], np.linalg.eigvalsh((D + D.conj().T) / 2)[0]

ok = True
lo, hi, dm = Tspec({(0,): 1.0}, 1)
print(f"(a) T_|0>: spectrum [{lo:+.6f}, {hi:+.6f}]"); ok &= abs(lo + 1 / 3) < 1e-12 and abs(hi - 1) < 1e-12
whi, wlo, wd = -9, 9, 9
for u in np.concatenate([np.linspace(0, 1, 21), [1.5, 3, 10]]):
    lo, hi, dm = Tspec({(1, 0): 1.0, (0, 1): u}, 2); whi = max(whi, hi); wlo = min(wlo, lo); wd = min(wd, dm)
print(f"(b) phi_u grid: max eig {whi:+.12f}, min eig {wlo:+.6f} (-11/9 = {-11/9:+.6f}), min lambda(D) {wd:.4f}")
ok &= whi < 1 + 1e-10 and wlo > -3 and wd > 0
rng = np.random.default_rng(7); cz = lambda *s: rng.normal(size=s) + 1j * rng.normal(size=s)
for k in range(3):
    x, y = cz(2), cz(2)
    lo, hi, dm = Tspec({(1, 0): x[0], (2, 0): x[1], (0, 1): y[0], (0, 2): y[1]}, 2)
    print(f"(c) random x, y: spectrum [{lo:+.6f}, {hi:+.6f}], min lambda(D) {dm:.4f}"); ok &= hi < 1 + 1e-10 and lo > -3 and dm > 0
worst = 9
for trial in range(8):
    n = [5, 6][trial % 2]
    U = [cz(2, 2), cz(2, 2)] + [np.outer(cz(2), [0, 1.0]) * rng.uniform(0.05, 2) for _ in range(n - 2)]
    if trial % 3 == 1: U[0] = 4 * U[0]
    if trial % 3 == 2: U[1] = np.outer(cz(2), cz(2))
    U = gauge(U)
    assert abs(np.linalg.eigvalsh(U[0].conj().T @ U[0] + U[1].conj().T @ U[1])[-1] - 1) < 1e-9
    for c in range(1, n + 1):
        worst = min(worst, np.linalg.eigvalsh(GeneralModel(U).ops(c)['K'])[0])
print(f"(d) A_0 + A_1 with eigenvalue 1, n = 5, 6, every c: min lambda_min(Kt_c) = {worst:+.3e}"); ok &= worst > -1e-12
print("PASS" if ok else "FAIL")
