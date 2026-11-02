"""Check the exact low thick critical points of R = q/||C||^2 on the rank-two variety (data/variety_critical_*.npy):
variety-criticality residual (P_A L and L P_B, L = Phi(C) - R C), R against its closed form, thickness, and the smallest
eigenvalue of the variety Hessian  q_r(A'B* + AB'*) + 2 Re<L, A'B'*>  (tangent metric, gauge directions removed).
usage: python exact_critical_points.py"""
import numpy as np
cases = [('data/variety_critical_n2_low.npy', 2, (9 - np.sqrt(41)) / 20, '(9 - sqrt 41)/20'),
         ('data/variety_critical_n3_golden.npy', 3, (3 - np.sqrt(5)) / 8, '(3 - sqrt 5)/8'),
         ('data/variety_critical_n3_low.npy', 3, (67 - np.sqrt(2329)) / 216, '(67 - sqrt 2329)/216')]
w = np.eye(3).reshape(9) / np.sqrt(3)
for fn, n, exact, label in cases:
    D = 3 ** n
    def Phi(C):
        X = C.reshape((3,) * (2 * n)).transpose([i for s in range(n) for i in (s, n + s)]).reshape((9,) * n)
        for s in range(n):
            tr = np.tensordot(w.conj(), X, axes=([0], [s])); X = X - 1.5 * np.moveaxis(np.tensordot(w, tr, axes=0), 0, s)
        return X.reshape([3, 3] * n).transpose([2 * s for s in range(n)] + [2 * s + 1 for s in range(n)]).reshape(D, D)
    def marg(P, s):
        T = P.reshape((3,) * (2 * n)); T = np.moveaxis(T, (s, n + s), (0, 1)).reshape(3, 3, 3 ** (n - 1), 3 ** (n - 1))
        return np.einsum('abxx->ab', T)
    C = np.load(fn); C = C / np.linalg.norm(C)
    r = np.vdot(C, Phi(C)).real; L = Phi(C) - r * C
    U, S, Vh = np.linalg.svd(C); Ua, Vb = U[:, :2], Vh[:2].conj().T
    res = np.linalg.norm(Ua.conj().T @ L) + np.linalg.norm(L @ Vb)
    thick = min(np.linalg.eigvalsh(marg(M @ M.conj().T, s))[0] for M in (Ua, Vb) for s in range(n))
    A = Ua * S[:2]; B = Vb
    # Hessian on the real basis of (A', B')
    T, meta = [], []
    for side in (0, 1):
        for k in range(D):
            for c in range(2):
                for ph in (1.0, 1j):
                    e = np.zeros(D); e[k] = 1
                    T.append(ph * np.outer(e, B[:, c].conj()) if side == 0 else np.conj(ph) * np.outer(A[:, c], e)); meta.append((side, k, c, ph))
    m = len(T); Tm = np.stack([t.reshape(-1) for t in T], 1); PT = np.stack([(Phi(t) - r * t).reshape(-1) for t in T], 1)
    H = (Tm.conj().T @ PT).real; G = (Tm.conj().T @ Tm).real
    half = m // 2
    for i in range(half):
        for j in range(half, m):
            _, ki, ci, pi = meta[i]; _, kj, cj, pj = meta[j]
            if ci == cj:
                v = (np.conj(L[ki, kj]) * pi * np.conj(pj)).real; H[i, j] += v; H[j, i] += v
    H = (H + H.T) / 2; G = (G + G.T) / 2
    wg, Ug = np.linalg.eigh(G); keep = wg > 1e-9 * wg.max(); P = Ug[:, keep] / np.sqrt(wg[keep])
    ev = np.linalg.eigvalsh(P.T @ H @ P)
    print(f"{fn}: n = {n}, R = {r:.12f} vs {label} = {exact:.12f} (diff {abs(r - exact):.1e}); variety-criticality residual "
          f"{res:.1e}; thickness {thick:.4f}; Hessian smallest {ev[0]:+.4f}, {np.sum(ev < -1e-8)} negative directions")
