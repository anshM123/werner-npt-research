"""Check of the two-copy point data/levi_positive_n2.npy: a thick rank-two C, critical for the two-sided local filters
C -> (g1 (x) g2) C (h1 (x) h2)^T, at which the Levi form q_r restricted to V_C = {XC + CY : X, Y 1-local} is positive
semidefinite (positive on the complement of C), while the real second variation of R on the orbit is still negative.
usage: python levi_positive_point.py [data/levi_positive_n2.npy]"""
import sys, numpy as np
n, t = 2, 0.5
D = 3 ** n
fn = sys.argv[1] if len(sys.argv) > 1 else 'data/levi_positive_n2.npy'
w = np.eye(3).reshape(9) / np.sqrt(3)
def Phi(C):
    X = C.reshape(3, 3, 3, 3).transpose(0, 2, 1, 3).reshape(9, 9)          # (a1 b1), (a2 b2)
    X = X - 3 * t * np.outer(w, w.conj() @ X)                              # copy 1
    X = X - 3 * t * np.outer(X @ w.conj(), w)                              # copy 2
    return X.reshape(3, 3, 3, 3).transpose(0, 2, 1, 3).reshape(D, D)
def q(C): return np.vdot(C, Phi(C)).real
def marg(M, s):
    T = M.reshape(3, 3, 3, 3)
    return np.einsum('ajbj->ab', T) if s == 0 else np.einsum('jajb->ab', T)
def lmul(E, C, s):
    T = C.reshape(3, 3, 3, 3)
    return (np.einsum('ab,bjcd->ajcd', E, T) if s == 0 else np.einsum('ab,jbcd->jacd', E, T)).reshape(D, D)
def rmul(C, E, s):
    T = C.reshape(3, 3, 3, 3)
    return (np.einsum('jkbd,ba->jkad', T, E) if s == 0 else np.einsum('jkcb,ba->jkca', T, E)).reshape(D, D)
C = np.load(fn); C = C / np.linalg.norm(C); r = q(C); L = Phi(C) - r * C
U, S, Vh = np.linalg.svd(C); A, B = U[:, :2], Vh[:2].conj().T
thick = min(np.linalg.eigvalsh(marg(Q @ Q.conj().T, s))[0] for Q in (A, B) for s in range(n))
crit = np.sqrt(sum(np.linalg.norm(marg(C @ L.conj().T, s)) ** 2 + np.linalg.norm(marg(L.conj().T @ C, s)) ** 2 for s in range(n)))
print(f"rank {np.sum(S > 1e-10 * S[0])}, s2/s1 = {S[1] / S[0]:.4f}, R = {r:.6f}, thickness = {thick:.4f}, two-sided gradient = {crit:.1e}")
# Levi form on V_C, with the C direction removed (rank-revealing)
dirs = []
for s in range(n):
    for a in range(3):
        for b in range(3):
            for ph in (1.0, 1j):
                E = np.zeros((3, 3), complex); E[a, b] = ph; dirs.append((s, E))
V = np.stack([lmul(E, C, s).reshape(-1) for s, E in dirs[::2]] + [rmul(C, E, s).reshape(-1) for s, E in dirs[::2]], 1)
c = C.reshape(-1); V = V - np.outer(c, c.conj() @ V)
Uv, Sv, _ = np.linalg.svd(V, full_matrices=False); Uv = Uv[:, Sv > 1e-9 * Sv[0]]
Qm = Uv.conj().T @ np.stack([Phi(Uv[:, j].reshape(D, D)).reshape(-1) for j in range(Uv.shape[1])], 1)
lev = np.linalg.eigvalsh((Qm + Qm.conj().T) / 2) - r
print(f"Levi form on the complement of C in V_C (dim {len(lev)}): smallest eigenvalues {np.round(lev[:4], 5)}")
# real second variation along e^{eX} C e^{eY}: q_r(XC + CY) + Re <L, X^2 C + 2 X C Y + C Y^2>  (X, Y 1-local, complex)
m = 2 * len(dirs); T = []; ops = []
for side in (0, 1):
    for s, E in dirs:
        T.append(lmul(E, C, s) if side == 0 else rmul(C, E, s)); ops.append((side, s, E))
PT = [Phi(Ti) - r * Ti for Ti in T]
def second(i, j):
    (si, ci, Ei), (sj, cj, Ej) = ops[i], ops[j]
    if si == 0 and sj == 0: return 0.5 * (lmul(Ei, lmul(Ej, C, cj), ci) + lmul(Ej, lmul(Ei, C, ci), cj))
    if si == 1 and sj == 1: return 0.5 * (rmul(rmul(C, Ei, ci), Ej, cj) + rmul(rmul(C, Ej, cj), Ei, ci))
    if si == 0: return rmul(lmul(Ei, C, ci), Ej, cj)
    return rmul(lmul(Ej, C, cj), Ei, ci)
H = np.zeros((m, m)); G = np.zeros((m, m))
for i in range(m):
    for j in range(i, m):
        H[i, j] = H[j, i] = np.vdot(T[i], PT[j]).real + np.vdot(L, second(i, j)).real
        G[i, j] = G[j, i] = np.vdot(T[i], T[j]).real
wg, Ug = np.linalg.eigh(G); keep = wg > 1e-8 * wg.max(); P = Ug[:, keep] / np.sqrt(wg[keep])
hev = np.linalg.eigvalsh(P.T @ H @ P)
print(f"real second variation of R on the two-sided orbit (tangent metric, dim {keep.sum()}): smallest {np.round(hev[:4], 4)}, "
      f"{np.sum(hev < -1e-8)} negative")
