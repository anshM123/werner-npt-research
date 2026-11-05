"""Exploration (floats): how much room is there in the splitting X = X' + Omega?
For each weight vector compute, on the soft space Omega = a w^T + w c^T (a perp w) [dim 2u - 1]:
  Fs   = soft form (exact F restricted to Omega),
  Cex  = exact Schur complement Fs - Bc^T F'^{-1} Bc (F' = F on X'-space),
  Cnorm= Fs - (12/5) Bc^T Pi' Bc  (hard part bounded by (5/12)||X'||^2),
  CL   = Fs - Bc^T (Lambda on V')^{-1} Bc (hard part bounded by the diagonal crossed-pair bound Lambda),
and report min generalized eigenvalue of (Cnorm, Fs), (CL, Fs), (Cex, Fs).  Positive = the bound closes."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cu_form import CU


def null_basis(A, tol=1e-12):
    u_, s_, vt = np.linalg.svd(A)
    r = (s_ > tol * max(1, s_.max())).sum()
    return vt[r:].T


def analyse(w, sig):
    u = len(w)
    cu = CU(list(w), sig, F=float)
    n = u * u
    # matrix of F on vec(X)
    Mf = np.zeros((n, n))
    basis = np.eye(n)
    def Fv(x):
        return cu.form(x.reshape(u, u).tolist())
    diag = [Fv(basis[i]) for i in range(n)]
    for i in range(n):
        Mf[i, i] = diag[i]
        for j in range(i + 1, n):
            Mf[i, j] = Mf[j, i] = (Fv(basis[i] + basis[j]) - diag[i] - diag[j]) / 2
    wv = np.asarray(w, float)
    Bw = null_basis(wv[None, :])                       # u x (u-1), orthonormal, perp w
    # X' basis: outer(Bw_i, Bw_j)
    Xp = np.array([np.outer(Bw[:, i], Bw[:, j]).ravel() for i in range(u - 1) for j in range(u - 1)]).T
    # soft basis: a w^T (a = Bw_i) and w c^T (c = e_j)
    Om = [np.outer(Bw[:, i], wv).ravel() for i in range(u - 1)] + [np.outer(wv, np.eye(u)[j]).ravel() for j in range(u)]
    Om = np.array(Om).T
    Fh = Xp.T @ Mf @ Xp
    Fs = Om.T @ Mf @ Om
    Bc = Xp.T @ Mf @ Om
    Cex = Fs - Bc.T @ np.linalg.solve(Fh, Bc)
    Cnorm = Fs - (12 / 5) * Bc.T @ Bc
    # Lambda bound
    s = wv ** 2
    d = np.array(cu.d)
    rho = s / d
    g = rho.sum()
    Lam = np.full((u, u), 0.5) - np.outer(rho, rho) / (2 + g)
    for h in range(u):
        Lam[h, h] = 1.5 - rho[h] + rho[h] ** 2 / (2 + g)
    Lm = Xp.T @ np.diag(Lam.ravel()) @ Xp
    CL = Fs - Bc.T @ np.linalg.solve(Lm, Bc)
    def gmin(C):
        L = np.linalg.cholesky(Fs)
        Li = np.linalg.inv(L)
        return np.linalg.eigvalsh(Li @ C @ Li.T)[0]
    return gmin(Cnorm), gmin(CL), gmin(Cex), np.linalg.eigvalsh(Fh)[0]


rng = np.random.default_rng(3)
for u in (2, 3, 4, 5, 6):
    worst = [np.inf] * 4
    for trial in range(150):
        spread = rng.choice([0.3, 1, 2, 4])
        w = np.exp(rng.uniform(-spread, spread, u))
        sig = np.exp(rng.uniform(-6, 6)) * (w ** 2).sum()
        r = analyse(w, sig)
        worst = [min(a, b) for a, b in zip(worst, r)]
    print(u, 'min gen.eig: norm-bound %.4f  Lambda-bound %.4f  exact %.4f ; min eig F on X-prime %.4f' % tuple(worst), flush=True)
