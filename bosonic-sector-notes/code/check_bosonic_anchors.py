"""Anchors inside the bosonic sector (NOTES.md, Section 6).
For p in Sym^n(C^3), M the form q_n on vec(C) (C in End(Sym^n)), V_p = p (x) I, D_p = V_p^* M V_p,
R_p = M V_p D_p^{-1} V_p^* M and Y_p = M - R_p.  Minimising q_n(|p><v| + |u><w|) over v in Sym^n leaves
(u (x) conj w)^* Y_p (u (x) conj w); p is a bosonic anchor if Y_p^Gamma >= 0 (Gamma: transpose on the row factor),
and then every bosonic X whose range lies in a plane containing p has q_n(X) >= 0, whatever its bosonic co-range.
Prints lambda_min(Y_p^Gamma) and the least value of Y_p on product vectors (alternating minimisation).
usage: python check_bosonic_anchors.py [n_max (default 4)]"""
import sys
import numpy as np
from bosonic import basis, q_and_grad

rng = np.random.default_rng(611)


def form_matrix(n):
    N = len(basis(n)); M = np.zeros((N * N, N * N), dtype=complex)
    for k in range(N * N):
        E = np.zeros((N, N), dtype=complex); E.flat[k] = 1
        _, G = q_and_grad(E, n)
        M[:, k] = G.reshape(-1)
    return (M + M.conj().T) / 2


def pt_row(Y, N):
    return Y.reshape(N, N, N, N).transpose(2, 1, 0, 3).reshape(N * N, N * N)


def stats(p, M, N):
    p = p / np.linalg.norm(p); V = np.kron(p.reshape(N, 1), np.eye(N))
    D = V.conj().T @ M @ V; Y = M - M @ V @ np.linalg.solve(D, V.conj().T @ M); Y = (Y + Y.conj().T) / 2
    lam = np.linalg.eigvalsh(pt_row(Y, N)).min()
    Yr = Y.reshape(N, N, N, N); best = np.inf
    for _ in range(6):
        u = rng.normal(size=N) + 1j * rng.normal(size=N)
        for _ in range(300):
            Au = np.einsum('i,ijkl,k->jl', u.conj(), Yr, u); e, U = np.linalg.eigh((Au + Au.conj().T) / 2); w = U[:, 0]
            Aw = np.einsum('j,ijkl,l->ik', w.conj(), Yr, w); e2, U2 = np.linalg.eigh((Aw + Aw.conj().T) / 2); u = U2[:, 0]
        best = min(best, e2[0])
    return lam, best


# sanity: the form matrix reproduces q on a random rank-two operator
n0 = 3; B = basis(n0); N0 = len(B); M0 = form_matrix(n0)
C = (rng.normal(size=(N0, 2)) + 1j * rng.normal(size=(N0, 2))) @ (rng.normal(size=(2, N0)) + 1j * rng.normal(size=(2, N0)))
print(f"form matrix check at n = 3: |vec(C)^* M vec(C) - q_n(C)| = {abs(np.vdot(C.reshape(-1), M0 @ C.reshape(-1)).real - q_and_grad(C, n0)[0]):.1e}")
for n in range(2, (int(sys.argv[1]) if len(sys.argv) > 1 else 4) + 1):
    B = basis(n); N = len(B); idx = {t: k for k, t in enumerate(B)}; M = form_matrix(n)
    tests = {'coherent |n,0,0>': [(n, 0, 0)], 'W |n-1,1,0>': [(n - 1, 1, 0)], '|n-2,1,1>': [(n - 2, 1, 1)],
             'Dicke |n/2,n/2,0>': [(n // 2, n - n // 2, 0)], 'GHZ': [(n, 0, 0), (0, n, 0), (0, 0, n)]}
    for name, occs in tests.items():
        p = np.zeros(N, complex)
        for o in occs: p[idx[o]] = 1
        lam, best = stats(p, M, N)
        print(f"n = {n}, {name:18s}: lambda_min(Y^Gamma) = {lam:+.3e}, least value on product vectors = {best:+.3e}", flush=True)
    for k in range(3):
        lam, best = stats(rng.normal(size=N) + 1j * rng.normal(size=N), M, N)
        print(f"n = {n}, random p {k}          : lambda_min(Y^Gamma) = {lam:+.3e}, least value on product vectors = {best:+.3e}", flush=True)
