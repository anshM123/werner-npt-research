"""Exploration (u >= 2): eliminate X_T exactly, remove the kernel (gauge w_U^T X_U = 0), eliminate the U-block
(X_U mod kernel, q).  The Schur complement G lives on the T-blocks (rho_t, zeta_t).  Coordinates per t:
(rho_t, z_t = w_U-component of zeta_t, z'_t = components perpendicular to w_U); rescale z_t by c = sqrt(1 + s_U/S_T).
Measure  A = max ||G_tt'||_2 S_T / sqrt(s_t s_t'),  B = min lambda_min(G_tt) / (1 - s_t/S_T),  and the weighted
dominance ratio  max_t sum_t' ||G_tt'|| sqrt(s_t'/s_t) / lambda_min(G_tt)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cf import TypeI


def null_basis(A, tol=1e-12):
    u, s, vt = np.linalg.svd(A)
    r = (s > tol * max(1, s.max())).sum()
    return vt[r:].T


def Gmat(w, T, rescale=True):
    ti = TypeI(w, T, F=float)
    M = np.array(ti.build())
    U = ti.U; u = len(U); w = np.asarray(w, float); s = w ** 2
    sT, sU = s[T].sum(), s[U].sum()
    wU = w[U]; nw = np.linalg.norm(wU)
    Bw = null_basis(wU[None, :])
    n = ti.n
    # change of basis: columns = new coordinates expressed in old ones
    cols = []; groups = {}
    # U-block: X_U mod kernel: X = X' (perp both sides) + a w^T (a perp w), plus q ; X_T free
    def XU(X):
        v = np.zeros(n)
        for a, r in enumerate(U):
            for b, h in enumerate(U):
                v[ti.idx[('P', r, h)]] = X[a, b]
        return v
    ublk = []
    for i in range(u - 1):
        for j in range(u - 1):
            ublk.append(XU(np.outer(Bw[:, i], Bw[:, j])))
    for i in range(u - 1):
        ublk.append(XU(np.outer(Bw[:, i], wU / nw)))
    v = np.zeros(n); v[ti.idx[('Q',)]] = 1; ublk.append(v)
    xt = []
    for t in T:
        for h in U:
            v = np.zeros(n); v[ti.idx[('P', t, h)]] = 1; xt.append(v)
    c = np.sqrt(1 + sU / sT) if rescale else 1.0
    tblk = {}
    for t in T:
        b = []
        v = np.zeros(n); v[ti.idx[('R', t)]] = 1; b.append(v)
        v = np.zeros(n)
        for j, h in enumerate(U):
            v[ti.idx[('Z', t, h)]] = c * wU[j] / nw
        b.append(v)
        for i in range(u - 1):
            v = np.zeros(n)
            for j, h in enumerate(U):
                v[ti.idx[('Z', t, h)]] = Bw[j, i]
            b.append(v)
        tblk[t] = b
    Bel = np.array(ublk + xt).T
    Tall = np.array([v for t in T for v in tblk[t]]).T
    A11 = Bel.T @ M @ Bel; A12 = Bel.T @ M @ Tall; A22 = Tall.T @ M @ Tall
    G = A22 - A12.T @ np.linalg.solve(A11, A12)
    G = (G + G.T) / 2
    d = u + 1
    return G, d, s, sT, sU, np.linalg.eigvalsh(A11)[0]


def measure(w, T):
    G, d, s, sT, sU, l11 = Gmat(w, T)
    m = len(T)
    lam = [np.linalg.eigvalsh(G[i * d:(i + 1) * d, i * d:(i + 1) * d])[0] for i in range(m)]
    A = 0; B = 9; R = 0
    for i, t in enumerate(T):
        B = min(B, lam[i] / (1 - s[t] / sT) if m > 1 else lam[i])
        acc = 0
        for j, tp in enumerate(T):
            if i == j: continue
            nr = np.linalg.norm(G[i * d:(i + 1) * d, j * d:(j + 1) * d], 2)
            A = max(A, nr * sT / np.sqrt(s[t] * s[tp]))
            acc += nr * np.sqrt(s[tp] / s[t])
        R = max(R, acc / lam[i])
    return A, B, R, np.linalg.eigvalsh(G)[0], l11


if __name__ == "__main__":
    rng = np.random.default_rng(11)
    for u in (2, 3):
        for m in (1, 2, 3, 4):
            k = u + m; T = list(range(m))
            Amax = 0; Bmin = 9; Rmax = 0; gmin = 9; lmin11 = 9
            for trial in range(150):
                w = np.exp(rng.uniform(-2.5, 2.5, k) * rng.choice([0.3, 1, 2]))
                if trial % 5 == 1: w[:m] *= np.exp(rng.uniform(-6, -1))
                if trial % 5 == 2: w[m:] *= np.exp(rng.uniform(-6, -1))
                A, B, R, g, l11 = measure(w, T)
                Amax = max(Amax, A); Bmin = min(Bmin, B); Rmax = max(Rmax, R); gmin = min(gmin, g); lmin11 = min(lmin11, l11)
            print(f"u={u} m={m}: sup A = {Amax:.4f}, inf B = {Bmin:.4f}, max dominance ratio = {Rmax:.4f}, min eig G = {gmin:.2e}, min eig eliminated block = {lmin11:.2e}", flush=True)
