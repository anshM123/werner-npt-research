"""Floating-point (numpy or mpmath) version of the closed forms of tlib.py, for exploration only.
Units 2^{k-1}, |s| = 1, copy t = 0 (sigma = x_0), U = 1..k-1.  States R_c (c in [k]), P_r (r in [k]), Q_c (c in U)."""
import numpy as np


def blocks(x):
    """x: squared weights, sum 1, copy 0 is t.  Returns dict of numpy blocks."""
    x = np.asarray(x, dtype=float)
    k = len(x)
    w = np.sqrt(x)
    sig = x[0]
    U = np.arange(1, k)
    N = np.array([1.0] + [2.0] * (k - 1))
    g = x / (2 * (1 - x))
    gam = g.sum()
    tau = (1 + g[0]) / (1 + gam / 2)
    xu = x[1:]
    gamp = (xu / (2 - sig - 2 * xu)).sum()
    taup = 1 / (1 + gamp / 2)
    DR = np.diag(N)
    taus = np.array([tau - 2] + [tau] * (k - 1))
    Lam = taus / (1 - x)
    P0 = np.eye(k) - np.outer(w, w)
    YRP = 0.5 * Lam[:, None] * P0
    Pi = np.diag(N) - np.outer(N * w, N * w) / (2 - sig)
    Dp = 1 / (2 - sig - 2 * xu)
    E = np.zeros((k, k - 1)); E[1:, :] = np.eye(k - 1)
    YPQ = taup * (np.outer(N * w, w[1:]) - (2 - sig) * E) * Dp[None, :]
    # N1
    dd = 1 - x[:, None] - x[None, :]
    N1 = np.outer(w, w) / (2 * np.where(np.eye(k) > 0, 1, dd))
    np.fill_diagonal(N1, 0)
    diag1 = np.array([sum(x[u] / (2 * (1 - x[c] - x[u])) for u in range(k) if u != c) for c in range(k)])
    N1 += np.diag(diag1)
    H = np.zeros((k, k - 1))
    for j, c in enumerate(U):
        H[c, j] = 1 + sig / (2 * (1 - sig - x[c]))
        H[0, j] = w[0] * w[c] / (2 * (1 - sig - x[c]))
    ddp = 2 - sig - 2 * xu[:, None] - 2 * xu[None, :]
    Np = np.outer(w[1:], w[1:]) / np.where(np.eye(k - 1) > 0, 1, ddp)
    np.fill_diagonal(Np, 0)
    diag2 = np.array([sum(xu[u] / (2 - sig - 2 * xu[c] - 2 * xu[u]) for u in range(k - 1) if u != c) for c in range(k - 1)])
    Np += np.diag(diag2)
    A1 = np.eye(k) + 0.5 * N1
    A2 = np.eye(k - 1) + 0.5 * Np
    YRQ = -np.linalg.solve(A1, H)
    YQQ = 2 * np.linalg.inv(A2)
    return dict(k=k, x=x, w=w, sig=sig, N=N, DR=DR, YRP=YRP, Pi=Pi, YPQ=YPQ, YRQ=YRQ, YQQ=YQQ, A1=A1, A2=A2, H=H,
                N1=N1, Np=Np, tau=tau, taup=taup, Lam=Lam, Dp=Dp, E=E, gam=gam, gamp=gamp)


def full(b):
    k = b['k']
    M = np.block([[b['DR'], b['YRP'], b['YRQ']], [b['YRP'].T, b['Pi'], b['YPQ']], [b['YRQ'].T, b['YPQ'].T, b['YQQ']]])
    return M


def phi(b):
    Bm = np.hstack([b['YRP'].T, b['YPQ']])        # P x (R, Q)
    C = np.block([[b['DR'], b['YRQ']], [b['YRQ'].T, b['YQQ']]])
    return C - Bm.T @ np.diag(1 / b['N']) @ Bm


def phip(b):
    D = np.block([[b['A1'], np.zeros((b['k'], b['k'] - 1))], [np.zeros((b['k'] - 1, b['k'])), b['A2']]])
    return D @ phi(b) @ D


def nmin(A):
    d = np.sqrt(np.diag(A))
    return np.linalg.eigvalsh(A / d[:, None] / d[None, :]).min()
