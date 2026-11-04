"""The refined concave-line step works where SC-infinity fails.

At the exact three-copy counterexample (W-type range, integer co-range; see check_sc_counterexamples.py), copy 0 has
no H with N_0(H) = 0 and Q_0(H) < 0.  The step of Theorem 3' only needs an indefinite H with F(H) < 0, where
F = Q_0 - R(C) N_0.  Along M = I + eps H, g(eps) = F(I + eps H) = 2 eps F(I, H) + eps^2 F(H) <= 2 eps F(I, H), so at the
end of the PSD interval where eps F(I, H) <= 0 we get q <= R(C) ||.||^2, with an end operator thin at copy 0.
Checks: F has negative directions with indefinite H; both ends are thin (local-support ratio ~ 1e-11); R at the
chosen end is below R(C)."""
import numpy as np
from scipy.linalg import sqrtm
from common import squeeze_forms, q, tangent_range, Bt

n = 3; D = 27
def ket(s):
    v = np.zeros(D); v[9 * s[0] + 3 * s[1] + s[2]] = 1; return v
u = tangent_range([np.eye(2)] * n, np.zeros(2), n).real
g1 = -ket((0,0,0)) - 3*ket((0,0,1)) - 3*ket((0,1,0)) + ket((0,2,0)) + 3*ket((1,0,0)) - ket((1,0,1)) - ket((1,1,2)) - ket((1,2,1)) + 2*ket((2,1,1))
g2 = ket((0,0,0)) - 3*ket((0,0,2)) - 3*ket((0,2,0)) - 2*ket((1,2,2)) + 3*ket((2,0,0)) + ket((2,0,2)) + ket((2,1,2)) + ket((2,2,1))
C = np.outer(u[0], g1) + np.outer(u[1], g2); C /= np.linalg.norm(C); R0 = q(C, n)
K, Kn = squeeze_forms(C, n, 0); nK = np.linalg.norm(Kn, 2)
lams = np.linspace(-6, 6, 12001); pd = [l for l in lams if np.linalg.eigvalsh(K - l * Kn)[0] > 0]
print(f"R(C) = {R0:.5f}; copy-0 pencil positive definite for lambda in about [{min(pd):.2f}, {max(pd):.2f}] (no Finsler direction)")
F = K - R0 * Kn; w, V = np.linalg.eigh(F)
H = np.einsum('k,kab->ab', V[:, 0], Bt); ev = np.linalg.eigvalsh(H)
FIH = np.einsum('k,kl,l->', np.array([np.trace(B.conj().T @ np.eye(3)).real for B in Bt]), F, V[:, 0])
print(f"F = Q_0 - R(C) N_0: smallest eigenvalue {w[0] / nK:+.4e} (relative to ||N_0||), eigenvector H with eigenvalues {np.round(ev, 3)}; F(I, H) = {FIH:+.3e}")
def squeeze(X, M):
    G = sqrtm(M); G = (G + G.conj().T) / 2; Gf = np.kron(G, np.eye(9)); return Gf @ X @ Gf
def ratio(X):
    T = np.moveaxis((X @ X.conj().T).reshape([3] * 6), [0, 3], [0, 3]).reshape(3, 9, 3, 9); e = np.linalg.eigvalsh(np.einsum('aibi->ab', T)); return e[0] / e[-1]
em, ep = -1 / ev[-1], -1 / ev[0]; res = []
for e in (em, ep):
    Y = squeeze(C, np.eye(3) + e * (1 - 1e-10) * H); res.append((e, q(Y / np.linalg.norm(Y), n), ratio(Y)))
    print(f"end eps = {e:+.4f}: R = {res[-1][1]:.5f}, local-support ratio at copy 0 = {res[-1][2]:.1e}")
star = res[0] if em * FIH <= 0 else res[1]
ok = w[0] < 0 and ev[0] < 0 < ev[-1] and star[1] < R0 and max(r[2] for r in res) < 1e-8
print(f"chosen end (eps F(I, H) <= 0): R = {star[1]:.5f} < R(C) = {R0:.5f}: {star[1] < R0}")
assert ok
print("PASS")
