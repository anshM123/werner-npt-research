"""Floating-point exploration of the map route of README.md, Section 5 (no claim rests on it).

P_a(X) = a Tr(X) I - (X + X^T)/4 on M_3, Choi matrix C_a = a I - (F + Omega Omega^*)/4.  For n copies:
  (1) m2(n) = min over Schmidt-rank-2 unit chi of <chi| C_a^{(x)n} |chi>, by alternating bottom eigenvectors with
      re-orthonormalisation of the two Schmidt pairs (random complex starts), compared with the product value
      (a - 3/4)(a - 1/2)^(n-1);
  (2) at a = 3/4 and a = 0.74: the same minimum (negative means P_a^{(x)n} is not 2-positive);
  (3) calibration: the map P_{9/10} (x) S, S with Choi matrix the rank-four projector I - sum of five product
      projectors (README, Section 6), whose tensor product with P_{9/10} is not 2-positive."""
import sys
import numpy as np

def choi_ops(a):
    e = np.eye(3)
    M = np.zeros((9, 9))
    for i in range(3):
        for j in range(3):
            Eij = np.outer(e[i], e[j]); M[:, 3 * i + j] = (a * np.trace(Eij) * np.eye(3) - (Eij + Eij.T) / 4).reshape(9)
    return M.reshape(3, 3, 3, 3)

def make_phi(site_maps):
    """site_maps: list of (3,3,3,3) superoperators (out_row, out_col, in_row, in_col), one per copy."""
    n = len(site_maps); D = 3 ** n
    def Phi(X):
        T = X.reshape([3] * (2 * n))
        for c, M4 in enumerate(site_maps):
            T = np.moveaxis(np.tensordot(M4, T, axes=([2, 3], [c, n + c])), [0, 1], [c, n + c])
        return T.reshape(D, D)
    def Phi_adj(Y):
        T = Y.reshape([3] * (2 * n))
        for c, M4 in enumerate(site_maps):
            T = np.moveaxis(np.tensordot(M4.conj(), T, axes=([0, 1], [c, n + c])), [0, 1], [c, n + c])
        return T.reshape(D, D)
    return Phi, Phi_adj, D

def margin(site_maps, starts, rng, iters=300):
    Phi, Phi_adj, D = make_phi(site_maps)
    def Lift(v, f):
        out = np.zeros((2 * D, 2 * D), dtype=complex)
        for r in range(2):
            for s in range(2):
                out[r * D:(r + 1) * D, s * D:(s + 1) * D] = f(np.outer(v[r], v[s].conj()))
        return (out + out.conj().T) / 2
    def orth(left, right):
        Ql, Rl = np.linalg.qr(left.T)
        U, sv, Wh = np.linalg.svd(Rl @ right, full_matrices=False)
        return (Ql @ U).T, np.diag(sv) @ Wh
    best = np.inf
    for st in range(starts):
        x = rng.normal(size=(2, D)) + 1j * rng.normal(size=(2, D))
        y = rng.normal(size=(2, D)) + 1j * rng.normal(size=(2, D))
        x, y = orth(x, y)
        val = np.inf
        for it in range(iters):
            w, V = np.linalg.eigh(Lift(x.conj(), Phi)); y = V[:, 0].reshape(2, D)
            y, x = orth(y, x)
            w2, V2 = np.linalg.eigh(Lift(y, Phi_adj)); x = V2[:, 0].reshape(2, D).conj()
            x, y = orth(x, y)
            if abs(w2[0] - val) < 1e-14:
                break
            val = w2[0]
        best = min(best, val)
    return best

rng = np.random.default_rng(2026)
if len(sys.argv) > 1 and sys.argv[1] == 'big':
    print('(4) five copies')
    for a in (0.8, 0.9):
        m = margin([choi_ops(a)] * 5, 3, rng)
        prod = (a - 0.75) * (a - 0.5) ** 4
        print(f'  a = {a}, n = 5: min found {m:.8f}, product value {prod:.8f}, ratio {m / prod:.6f}')
    m = margin([choi_ops(0.75)] * 5, 3, rng)
    print(f'  a = 0.75, n = 5: min over Schmidt-rank-2 = {m:.3e}')
    sys.exit(0)
print('(1) multiplicative margin')
for a in (0.8, 0.9):
    for n in (1, 2, 3, 4):
        m = margin([choi_ops(a)] * n, 14 if n < 4 else 4, rng)
        prod = (a - 0.75) * (a - 0.5) ** (n - 1)
        print(f'  a = {a}, n = {n}: min found {m:.8f}, product value {prod:.8f}, ratio {m / prod:.6f}')
print('(2) the threshold')
for a, ns in ((0.75, (2, 3, 4)), (0.74, (1, 2))):
    for n in ns:
        m = margin([choi_ops(a)] * n, 6 if n < 4 else 3, rng)
        print(f'  a = {a}, n = {n}: min over Schmidt-rank-2 = {m:.3e}')
print('(3) calibration on P_{9/10} (x) S')
e = np.eye(3)
v = [np.kron(e[0], e[0] - e[1]), np.kron(e[0] - e[1], e[2]), np.kron(e[2], e[1] - e[2]), np.kron(e[1] - e[2], e[0]),
     np.kron(e[0] + e[1] + e[2], e[0] + e[1] + e[2])]
rho = np.eye(9) - sum(np.outer(x, x) / (x @ x) for x in v)
Sm = np.zeros((3, 3, 3, 3))
C4 = rho.reshape(3, 3, 3, 3)                       # <i k| C |j l>
for i in range(3):
    for j in range(3):
        Eij = np.outer(e[i], e[j]); Sm[:, :, i, j] = np.einsum('ji,ikjl->kl', Eij, C4)
m = margin([choi_ops(0.9), Sm], 10, rng, iters=2000)
print(f'  min over Schmidt-rank-2 = {m:.3e} (negative: the search finds this known failure)')
