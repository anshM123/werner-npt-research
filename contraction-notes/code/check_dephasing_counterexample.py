"""The dephasing condition fails at two copies (exact certificate).

C = Y Z^* on C^3 (x) C^3 with Gaussian-integer factors. For a copy s and an orthonormal basis (u_i) of it, the
dephasing defect is  sum_i q_1(<u_i|_s C |u_i>_s) - q_1(Tr_s C)/2.  If H_1, H_2 are Hilbert-Schmidt orthonormal and
span the traceless diagonal matrices of the basis, then  defect = sum_j q_1(Tr_s[(H_j (x) I) C]) - q_1(Tr_s C)/6.
Writing H_j = sum_a r_{j,a} K_a in the Hermitian traceless basis K_a below, the vectors r_j satisfy
sum_j r_j^T G r_j = 2 and sum_j r_j^T D_c r_j = 0 (from H_1^2 + H_2^2 = (2/3) I), so for every basis
    defect <= 2 lambda - beta_s   whenever   S = lambda G - A_s + sum_c y_c D_c  is positive semidefinite.
Here A_s[a, b] = Re q_1(Z_a, Z_b) with Z_a = Tr_s[(K_a (x) I) C], beta_s = q_1(Tr_s C)/6, D_c[a, b] = Re Tr K_c(K_a K_b +
K_b K_a)/2. Positivity of S is certified by diagonal dominance of the integer matrix 2S. Everything is exact (SymPy).
A direct numerical maximisation of the defect over bases is printed as a cross-check."""
import numpy as np
import sympy as sp
from scipy.linalg import expm
from scipy.optimize import minimize

I = sp.I
Y = sp.Matrix([[-I, -1 + 4*I], [2 - 2*I, 1 - 4*I], [-1 - I, -1 - I], [3*I, -4*I], [0, -1 + 3*I], [-2*I, -1 - I],
               [-3 + 2*I, 1 - I], [1 + 2*I, I], [-1, I]])
Z = sp.Matrix([[1 + 3*I, 4 + 5*I], [-1 + I, -1 + I], [2 + 2*I, 1], [1 + 7*I, -5 - 6*I], [3 + 4*I, 2 + 3*I],
               [2 + 2*I, 2 + I], [-8 + 6*I, 1 - 4*I], [-3 + 3*I, -3 + 2*I], [4 - 4*I, 1 - 2*I]])
C = Y * Z.H                                   # rows/columns ordered 00, 01, ..., 22


def hs(A, B):
    return sp.expand(sum(sp.conjugate(A[i, j]) * B[i, j] for i in range(A.rows) for j in range(A.cols)))


def ptr(M, s):
    """Partial trace over copy s of a 9 x 9 operator."""
    if s == 0:
        return sp.Matrix(3, 3, lambda a, b: sum(M[3*k + a, 3*k + b] for k in range(3)))
    return sp.Matrix(3, 3, lambda a, b: sum(M[3*a + k, 3*b + k] for k in range(3)))


def q1(A, B=None):
    B = A if B is None else B
    return sp.expand(hs(A, B) - sp.Rational(1, 2) * sp.conjugate(A.trace()) * B.trace())


T0, T1 = ptr(C, 0), ptr(C, 1)
norm2 = hs(C, C)
q2 = sp.nsimplify(sp.expand(norm2 - (hs(T0, T0) + hs(T1, T1)) / 2 + sp.Abs(C.trace())**2 / 4))
print(f"rank Y = {Y.rank()}, rank Z = {Z.rank()};  ||C||^2 = {norm2};  q_2(C) = {q2}  (ratio {float(q2 / norm2):.8f})")

E = lambda i, j: sp.Matrix(3, 3, lambda a, b: 1 if (a, b) == (i, j) else 0)
K = [E(0, 1) + E(1, 0), I*(E(1, 0) - E(0, 1)), E(0, 2) + E(2, 0), I*(E(2, 0) - E(0, 2)),
     E(1, 2) + E(2, 1), I*(E(2, 1) - E(1, 2)), sp.diag(1, -1, 0), sp.diag(1, 1, -2)]
G = sp.Matrix(8, 8, lambda a, b: sp.re(sp.expand((K[a].H * K[b]).trace())))
assert G == sp.diag(2, 2, 2, 2, 2, 2, 2, 6)
D = [sp.Matrix(8, 8, lambda a, b: sp.re(sp.expand((K[c] * (K[a]*K[b] + K[b]*K[a]) / 2).trace()))) for c in range(8)]
certificate = {0: (1267, [396, -100, 56, -180, -328, -708, 124, -100]),
               1: (1276, [152, 196, 172, -56, 180, 160, -180, -408])}
for s in (0, 1):
    Kon = lambda Kc: sp.kronecker_product(Kc, sp.eye(3)) if s == 0 else sp.kronecker_product(sp.eye(3), Kc)
    Za = [ptr(Kon(K[a]) * C, s) for a in range(8)]
    A = sp.Matrix(8, 8, lambda a, b: sp.re(q1(Za[a], Za[b])))
    beta = sp.re(q1(ptr(C, s))) / 6
    lam, y = certificate[s]
    S2 = 2 * (lam * G - A + sum((y[c] * D[c] for c in range(8)), sp.zeros(8, 8)))
    assert all(S2[a, b].is_integer for a in range(8) for b in range(8)) and S2 == S2.T
    surplus = [S2[a, a] - sum(abs(S2[a, b]) for b in range(8) if b != a) for a in range(8)]
    assert min(surplus) >= 0
    print(f"copy {s}: beta = {beta};  2S is an integer matrix, diagonal-dominance surpluses {surplus};"
          f"  so the defect is <= 2*{lam} - beta = {2*lam - beta} in every orthonormal basis")

# cross-check: maximise the defect numerically over bases (exp of a Hermitian matrix)
Cn = np.array(C, dtype=complex)


def ptr_n(M, s):
    T = M.reshape(3, 3, 3, 3)
    return np.einsum('kakb->ab', T) if s == 0 else np.einsum('akbk->ab', T)


def q1n(X):
    return np.vdot(X, X).real - 0.5 * abs(np.trace(X))**2


def defect(U, s):
    tot = 0.0
    for i in range(3):
        P = np.outer(U[:, i], U[:, i].conj())
        tot += q1n(ptr_n((np.kron(P, np.eye(3)) if s == 0 else np.kron(np.eye(3), P)) @ Cn, s))
    return tot - q1n(ptr_n(Cn, s)) / 2


def unitary(h):
    H = np.zeros((3, 3), complex)
    iu = np.triu_indices(3, 1)
    H[np.diag_indices(3)] = h[:3]
    H[iu] = h[3:6] + 1j * h[6:9]
    return expm(1j * (H + np.triu(H, 1).conj().T))


rng = np.random.default_rng(0)
for s in (0, 1):
    best = max(-minimize(lambda h: -defect(unitary(h), s), rng.normal(size=9) * 2, method='BFGS').fun
               for _ in range(100))
    print(f"copy {s}: numerical maximum of the defect over 100 started bases = {best:.4f}")
