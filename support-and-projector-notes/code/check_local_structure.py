"""Checks for Section 6 of NOTES.md (local structure at zero planes), numerically at small n:
V1 projector-core zeros Z = P (x) ab^*: q(Z) = 0, L(Z) = -|2><2| (x) L_k(ab^*), and q(C(t)) = t^2 q_k(xb^* + ay^*) (k = 1, 2);
V2 spin-flip wedge zeros C0 = a(Theta b)^* - b(Theta a)^*: q = 0, stationarity, nonnegative Hessian on rank-two matrices,
   kernel dimension 4*2^n - 6 + 4n (+ 8 gauge) at the GHZ zero and at a generic pair (n = 1, 2);
V4 m(L) <= 2 sigma(L) on random planes and planes near the GHZ plane (n = 2)."""
import numpy as np, itertools
rng = np.random.default_rng(7)
def M1():
    Om = np.zeros(9); Om[[0, 4, 8]] = 1; return np.eye(9) - 0.5 * np.outer(Om, Om)
def Mn(n):
    m = M1(); out = np.array([[1.0]])
    for _ in range(n): out = np.kron(out, m)
    return out
def vecC(C, n):
    """vec in grouped (row_1 col_1)(row_2 col_2)... order"""
    T = C.reshape([3] * n + [3] * n)
    perm = [x for s in range(n) for x in (s, n + s)]
    return T.transpose(perm).reshape(-1)
def unvec(v, n):
    T = v.reshape([3] * (2 * n)); inv = list(range(0, 2 * n, 2)) + list(range(1, 2 * n, 2))
    return T.transpose(inv).reshape(3 ** n, 3 ** n)
def q(C, n, M): v = vecC(C, n); return np.real(np.conj(v) @ M @ v)
def Lop(C, n, M): return unvec(M @ vecC(C, n), n)
def basis(n, s):
    v = np.zeros(3 ** n); i = 0
    for c in s: i = 3 * i + c
    v[i] = 1; return v
# ---------- V1: projector-core zeros ----------
for k in (1, 2):
    n = k + 1; M = Mn(n); Mk = Mn(k)
    a = rng.normal(size=3 ** k) + 1j * rng.normal(size=3 ** k); b = rng.normal(size=3 ** k) + 1j * rng.normal(size=3 ** k)
    P = np.diag([1, 1, 0]).astype(complex)
    Z = np.kron(P, np.outer(a, b.conj()))
    LZ = Lop(Z, n, M); pred = -np.kron(np.diag([0, 0, 1]), Lop(np.outer(a, b.conj()), k, Mk))
    x = rng.normal(size=3 ** k) + 1j * rng.normal(size=3 ** k); y = rng.normal(size=3 ** k) + 1j * rng.normal(size=3 ** k)
    E01 = np.zeros((3, 3)); E01[0, 1] = 1
    t = 0.37; Ct = Z + t * np.kron(E01, np.outer(x, b.conj()) + np.outer(a, y.conj()))
    print(f"V1 k={k}: q(Z) = {q(Z, n, M):.1e}; |L(Z) - pred| = {np.abs(LZ - pred).max():.1e}; rank C(t) = {np.linalg.matrix_rank(Ct, 1e-9)}; "
          f"q(C(t)) - t^2 q_k(xb*+ay*) = {q(Ct, n, M) - t**2 * q(np.outer(x, b.conj()) + np.outer(a, y.conj()), k, Mk):.1e}")
# ---------- V2: spin-flip wedge zeros ----------
def theta(v, n):
    """Theta = J^{(x)n} complex conjugation on Q^{(x)n}, J|0> = |1>, J|1> = -|0>; v given on (C^3)^{(x)n} supported on {0,1}^n"""
    J = np.zeros((3, 3)); J[1, 0] = 1; J[0, 1] = -1
    Jn = np.array([[1.0]])
    for _ in range(n): Jn = np.kron(Jn, J)
    return Jn @ v.conj()
def qubitvec(n, rng, kind='rand'):
    v = np.zeros(3 ** n, complex)
    for s in itertools.product([0, 1], repeat=n):
        i = 0
        for c in s: i = 3 * i + c
        v[i] = rng.normal() + 1j * rng.normal()
    return v
def hessian(X0, Y0, n, M):
    """real quadratic form H(A,B) = q(A Y0^* + X0 B^*) + 2 Re<L(C0), A B^*>, (A,B) in C^{N x 2} x C^{N x 2} as real vector"""
    N = 3 ** n; C0 = X0 @ Y0.conj().T; G = Lop(C0, n, M)
    dim = 8 * N
    def unpack(z):
        A = (z[:2 * N] + 1j * z[2 * N:4 * N]).reshape(N, 2); B = (z[4 * N:6 * N] + 1j * z[6 * N:]).reshape(N, 2); return A, B
    # linear part: delta(z) as complex matrix -> realified columns
    cols = []
    for i in range(dim):
        z = np.zeros(dim); z[i] = 1; A, B = unpack(z)
        cols.append(vecC(A @ Y0.conj().T + X0 @ B.conj().T, n))
    Lm = np.array(cols).T                      # complex 9^n x dim
    H = np.real(Lm.conj().T @ M @ Lm)
    # bilinear curvature part: 2 Re <G, A B^*> = 2 Re Tr(G^* A B^*)
    Cm = np.zeros((dim, dim))
    for i in range(4 * N):
        zi = np.zeros(dim); zi[i] = 1; A, _ = unpack(zi)
        for j in range(4 * N, dim):
            zj = np.zeros(dim); zj[j] = 1; _, B = unpack(zj)
            val = 2 * np.real(np.trace(G.conj().T @ A @ B.conj().T))
            Cm[i, j] += val / 2; Cm[j, i] += val / 2
    return (H + H.T) / 2 + Cm
for n in (1, 2):
    M = Mn(n); N = 3 ** n
    for kind in ('generic', 'GHZ'):
        if kind == 'GHZ':
            a = basis(n, (0,) * n).astype(complex); b = basis(n, (1,) * n).astype(complex)
        else:
            a = qubitvec(n, rng); b = qubitvec(n, rng)
        X0 = np.column_stack([a, b]); Y0 = np.column_stack([theta(b, n), -theta(a, n)])
        C0 = X0 @ Y0.conj().T
        G = Lop(C0, n, M)
        # stationarity: G orthogonal to tangent {A Y0^* + X0 B^*}
        Ur, _, _ = np.linalg.svd(X0, full_matrices=False); Vr, _, _ = np.linalg.svd(Y0, full_matrices=False)
        Pr = Ur @ Ur.conj().T; Pc = Vr @ Vr.conj().T
        tangent_part = G - (np.eye(N) - Pr) @ G @ (np.eye(N) - Pc)
        H = hessian(X0, Y0, n, M); ev = np.linalg.eigvalsh(H)
        ker = int(np.sum(np.abs(ev) < 1e-9)); d = 2 ** n
        print(f"V2 n={n} {kind}: q(C0) = {q(C0, n, M):.1e}; |tangent part of L(C0)| = {np.abs(tangent_part).max():.1e}; "
              f"Hessian min eig {ev[0]:+.2e}; kernel dim {ker} (GHZ prediction {4*d - 6 + 4*n + 8})")
# ---------- V4: m(L) <= 2 sigma(L) ----------
def margin_and_sigma(U, n, M):
    N = 3 ** n
    Q = np.zeros((2 * N, 2 * N), complex)
    for j in range(2):
        for xc in range(N):
            Cm = np.zeros((N, N), complex); Cm[:, xc] = U[:, j]
            col = vecC(Cm, n); MC = unvec(M @ col, n)
            Q[:, j * N + xc] = (U.conj().T @ MC).reshape(-1)
    m = np.linalg.eigvalsh((Q + Q.conj().T) / 2)[0]
    Pi = U @ U.conj().T; T = Pi.reshape([3] * (2 * n)); sig = 0; tau = 0
    for v in range(n):
        axes = [i for i in range(n) if i != v]
        rho = np.trace(T, axis1=axes[0], axis2=n + axes[0]) if n == 2 else None
        if n == 2:
            rho = np.einsum('aibi->ab' if v == 0 else 'iaib->ab', T.reshape(3, 3, 3, 3))
        lm = np.linalg.eigvalsh(rho)[0]; sig += lm; tau = max(tau, lm)
    return m, sig, tau
n = 2; M = Mn(n); worst = 0
for trial in range(300):
    if trial % 2 == 0:
        U, _ = np.linalg.qr(rng.normal(size=(9, 2)) + 1j * rng.normal(size=(9, 2)))
    else:
        g = np.zeros((9, 2), complex); g[0, 0] = 1; g[4, 1] = 1     # GHZ plane span{|00>,|11>}
        U, _ = np.linalg.qr(g + 10 ** rng.uniform(-3, -0.5) * (rng.normal(size=(9, 2)) + 1j * rng.normal(size=(9, 2))))
    m, sig, tau = margin_and_sigma(U, n, M)
    worst = max(worst, m / (2 * sig))
print(f"V4 n=2: max m(L)/(2 sigma(L)) over 300 planes (random and near GHZ) = {worst:.4f} (claim <= 1)")
