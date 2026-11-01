"""Numerical checks of claims 1, 2, 4 and of the identities behind claim 3 (binary merger, normal form).
usage: python verify_claims.py   (NumPy only, seconds)"""
import numpy as np
from bm_core import G, pt, RHO, OMEGA, D, N, rand_c

g = np.random.default_rng(2026)
e = np.eye(D)

# ---- claim 1: A = |a><i|, B = |b><j|, i != j  ->  G = |ab><ab|
err1 = 0.0
for i in range(D):
    for j in range(D):
        a, b = rand_c(g, D), rand_c(g, D)
        A, B = np.outer(a, e[i]), np.outer(b, e[j])
        v = np.kron(a, b)
        ref = np.outer(v, v.conj()) * (0.5 if i == j else 1.0)
        err1 = max(err1, np.abs(G(A, B) - ref).max() / np.abs(ref).max())
print(f"claim 1: G(|a><i|, |b><j|) = |ab><ab| (i != j), = |ab><ab|/2 (i = j); max rel. error {err1:.1e}")

# ---- claim 2: G(A,B)^Gamma = AA^* (x) conj(BB^*) - (1/2) v v^*,  v = (A (x) conj B)|Omega> = vec(A B^*)
err2 = 0.0; err2b = 0.0
for _ in range(200):
    A, B = rand_c(g, D, D), rand_c(g, D, D)
    v = np.kron(A, B.conj()) @ OMEGA
    ref = np.kron(A @ A.conj().T, (B @ B.conj().T).conj()) - 0.5 * np.outer(v, v.conj())
    err2 = max(err2, np.abs(pt(G(A, B)) - ref).max() / np.abs(ref).max())
    err2b = max(err2b, np.abs(v - (A @ B.conj().T).reshape(N)).max())
print(f"claim 2: PT formula max rel. error {err2:.1e}; v = vec(A B^*) error {err2b:.1e}")
# number of negative PT eigenvalues of an atom with invertible filters, and the 'Lorentz' identity <v|(P(x)Q^T)^{-1}|v> = 3
negs = []; lor = []
for _ in range(200):
    A, B = rand_c(g, D, D), rand_c(g, D, D)
    ev = np.linalg.eigvalsh(pt(G(A, B)))
    negs.append(int(np.sum(ev < -1e-12 * abs(ev).max())))
    v = (A @ B.conj().T).reshape(N)
    PQ = np.kron(A @ A.conj().T, (B @ B.conj().T).conj())
    lor.append(np.vdot(v, np.linalg.solve(PQ, v)).real)
print(f"   atoms with random filters: #negative PT eigenvalues = {sorted(set(negs))}; <v|(P (x) conj Q)^-1|v> in "
      f"[{min(lor):.12f}, {max(lor):.12f}] (= 3: every atom has the same 'negativity strength')")

# ---- claim 4: Tr(W G(A,B)) = p_W(A,B) = vec(A)^* M_W(B) vec(A)
R4 = RHO.reshape(D, D, D, D)
def M_W(W, B):
    W4 = W.reshape(D, D, D, D)
    M = np.einsum('ijkl,Yl,XYxy,yj->Xkxi', R4, B.conj(), W4, B)
    return M.reshape(N, N)
def p_W(W, A, B):
    a = [A[:, i] for i in range(D)]; b = [B[:, j] for j in range(D)]
    s = 0
    for i in range(D):
        for j in range(D):
            s += np.vdot(np.kron(a[i], b[j]), W @ np.kron(a[i], b[j]))
            s -= 0.5 * np.vdot(np.kron(a[i], b[j]), W @ np.kron(a[j], b[i]))
    return s
err4 = 0.0; err4b = 0.0; herm = 0.0
for _ in range(200):
    W = rand_c(g, N, N); W = W + W.conj().T
    A, B = rand_c(g, D, D), rand_c(g, D, D)
    t = np.trace(W @ G(A, B))
    err4 = max(err4, abs(t - p_W(W, A, B)) / abs(t))
    M = M_W(W, B)
    err4b = max(err4b, abs(t - np.vdot(A.reshape(N), M @ A.reshape(N))) / abs(t))
    herm = max(herm, np.abs(M - M.conj().T).max())
    # quadratic (degree (1,1)) in B
    M2 = M_W(W, 2.0 * B)
    err4b = max(err4b, np.abs(M2 - 4 * M).max() / np.abs(M).max())
print(f"claim 4: Tr(W G) = p_W(A,B): max rel. error {err4:.1e}; = vec(A)^* M_W(B) vec(A) with M_W Hermitian "
      f"(asym {herm:.1e}) and quadratic in B: {err4b:.1e}")
# M_W(B) for W = I and W = rho^Gamma-type witnesses: W = I is in the dual cone (M_I(B) >= 0); the flip F is not.
mins_I = min(np.linalg.eigvalsh(M_W(np.eye(N), rand_c(g, D, D)))[0] for _ in range(200))
from bm_core import F9
mins_F = min(np.linalg.eigvalsh(M_W(F9, rand_c(g, D, D)))[0] for _ in range(200))
print(f"   min eigenvalue of M_I(B) over 200 random B: {mins_I:+.3e} (>= 0);  of M_F(B): {mins_F:+.3e} (F not in K*)")

# ---- claim 3 identities
# (i) normal form: G(A1,B1) + G(A2,B2) = (A1 (x) B1) [rho + G(A1^-1 A2, B1^-1 B2)] (A1 (x) B1)^*
err3 = 0.0
for _ in range(100):
    A1, B1, A2, B2 = (rand_c(g, D, D) for _ in range(4))
    K = np.kron(A1, B1)
    lhs = G(A1, B1) + G(A2, B2)
    rhs = K @ (RHO + G(np.linalg.solve(A1, A2), np.linalg.solve(B1, B2))) @ K.conj().T
    err3 = max(err3, np.abs(lhs - rhs).max() / np.abs(lhs).max())
print(f"claim 3 (normal form identity): max rel. error {err3:.1e}")
# (ii) singular filters give separable atoms: rank A <= 2 -> G(A,B) is (A(x)B) sigma (A(x)B)^* with sigma 2 x 3 PPT
worst = np.inf
for _ in range(200):
    A = rand_c(g, D, 2) @ rand_c(g, 2, D); B = rand_c(g, D, D)
    Y = G(A, B); worst = min(worst, np.linalg.eigvalsh(pt(Y))[0] / np.trace(Y).real)
print(f"   rank-2 filter atoms: min PT eigenvalue / Tr = {worst:+.2e} (PPT; separable via the 2 x 3 compression of rho)")
# (iii) symmetries of the normal form: rho + G(C, D) -> (U(x)U)(.)(U(x)U)^* = rho + G(UC, UD); G(CV, DV) = G(C, D)
C, Dm = rand_c(g, D, D), rand_c(g, D, D)
U = np.linalg.qr(rand_c(g, D, D))[0]; V = np.linalg.qr(rand_c(g, D, D))[0]
KU = np.kron(U, U)
e1 = np.abs(KU @ (RHO + G(C, Dm)) @ KU.conj().T - (RHO + G(U @ C, U @ Dm))).max()
e2 = np.abs(G(C @ V, Dm @ V) - G(C, Dm)).max()
print(f"   symmetries: U(x)U covariance {e1:.1e}; right unitary invariance G(CV, DV) = G(C, D): {e2:.1e}")
# (iv) PT quadratic form: <z|G(A,B)^Gamma|z> = phi(A^* Z B), phi(Y) = ||Y||^2 - |Tr Y|^2 / 2, z = vec(Z)
err5 = 0.0
for _ in range(100):
    A, B, Z = rand_c(g, D, D), rand_c(g, D, D), rand_c(g, D, D)
    z = Z.reshape(N); Y = A.conj().T @ Z @ B
    val = np.vdot(z, pt(G(A, B)) @ z).real
    err5 = max(err5, abs(val - (np.linalg.norm(Y) ** 2 - abs(np.trace(Y)) ** 2 / 2)) / abs(val))
print(f"   <z|G(A,B)^Gamma|z> = ||A^* Z B||^2 - |Tr(A^* Z B)|^2 / 2: max rel. error {err5:.1e}")
