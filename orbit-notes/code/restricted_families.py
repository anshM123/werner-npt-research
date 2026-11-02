"""Checks of three-copy squeeze-critical points at which a natural family of squeeze directions does not descend.
For C (rank two, thick) and G(x) = (x)_s exp(H_s/2), H_s = sum_k x_{s,k} E_k (E_k an orthonormal basis of Herm(3)),
f(x) = R_{1/2}(G C G).  The script prints R, the thickness, |grad f(0)| and the smallest eigenvalue of the Hessian of f
at 0 restricted to:
  torus   - traceless squeezes diagonal in the eigenbasis of mu_s = Tr_other(CC^* + C^*C) at every copy (2n directions);
  su(K)   - squeezes in the plane K_s of the two smaller mu_s-eigenvectors, su(2) of K_s at every copy (3n directions);
  pairs   - all squeezes at two copies (for each pair; identity directions removed);
  full    - all squeezes (identity directions removed).
Derivatives by central differences (step 1e-4; error ~1e-8).
usage: python restricted_families.py data/torus_pd_n3.npy [more .npy files]"""
import sys, itertools, numpy as np
from scipy.linalg import expm
t = 0.5
basis = []
for i in range(3):
    for j in range(3):
        M = np.zeros((3, 3), complex)
        if i == j: M[i, i] = 1
        elif i < j: M[i, j] = M[j, i] = 1 / np.sqrt(2)
        else: M[i, j] = 1j / np.sqrt(2); M[j, i] = -1j / np.sqrt(2)
        basis.append(M)
basis = np.array(basis)
def setup(n):
    D = 3 ** n; w = np.eye(3).reshape(9) / np.sqrt(3)
    def Phi(C):
        X = C.reshape((3,) * (2 * n)).transpose([i for s in range(n) for i in (s, n + s)]).reshape((9,) * n)
        for s in range(n):
            tr = np.tensordot(w.conj(), X, axes=([0], [s])); X = X - 3 * t * np.moveaxis(np.multiply.outer(w, tr), 0, s)
        return X.reshape([3, 3] * n).transpose([2 * s for s in range(n)] + [2 * s + 1 for s in range(n)]).reshape(D, D)
    return D, Phi
def kron_all(ms):
    out = np.ones((1, 1))
    for M in ms: out = np.kron(out, M)
    return out
def marg(M, s, n):
    T = M.reshape((3,) * (2 * n)); T = np.moveaxis(T, (s, n + s), (0, 1)).reshape(3, 3, 3 ** (n - 1), 3 ** (n - 1))
    return np.einsum('abxx->ab', T)
for fn in sys.argv[1:]:
    C = np.load(fn); n = int(round(np.log(C.shape[0]) / np.log(3))); D, Phi = setup(n)
    C = C / np.linalg.norm(C)
    R = lambda X: np.vdot(X, Phi(X)).real / np.vdot(X, X).real
    def f(x):
        G = kron_all([expm(0.5 * np.einsum('k,kab->ab', x[9 * s:9 * s + 9], basis)) for s in range(n)])
        return R(G @ C @ G)
    N = 9 * n; h = 1e-4; f0 = f(np.zeros(N))
    grad = np.array([(f(h * e) - f(-h * e)) / (2 * h) for e in np.eye(N)])
    H = np.zeros((N, N))
    for i in range(N):
        for j in range(i, N):
            ei, ej = np.eye(N)[i] * h, np.eye(N)[j] * h
            H[i, j] = H[j, i] = (f(ei + ej) - f(ei - ej) - f(-ei + ej) + f(-ei - ej)) / (4 * h * h)
    U, S, Vh = np.linalg.svd(C); A, B = U[:, :2], Vh[:2].conj().T
    thick = min(np.linalg.eigvalsh(marg(Q @ Q.conj().T, s, n))[0] for Q in (A, B) for s in range(n))
    Mu = C @ C.conj().T + C.conj().T @ C; Vs = [np.linalg.eigh(marg(Mu, s, n))[1] for s in range(n)]
    coords = lambda s, Hs: np.concatenate([np.einsum('kab,ba->k', basis, Hs).real if u == s else np.zeros(9) for u in range(n)])
    def restricted(cols):
        Q, _ = np.linalg.qr(np.array(cols).T); return np.linalg.eigvalsh(Q.T @ H @ Q)[0]
    torus = [coords(s, Vs[s] @ np.diag(v) @ Vs[s].conj().T) for s in range(n) for v in ([1, -1, 0], [1, 1, -2])]
    pauli = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])]
    suk = [coords(s, Vs[s][:, :2] @ p @ Vs[s][:, :2].conj().T) for s in range(n) for p in pauli]
    tl = lambda s: [coords(s, E) for E in basis if abs(np.trace(E)) < 1e-12] + [coords(s, np.diag(v)) for v in ([1, -1, 0], [1, 1, -2])]
    pairs = [restricted(tl(a) + tl(b)) for a, b in itertools.combinations(range(n), 2)]
    full = restricted(sum((tl(s) for s in range(n)), []))
    print(f"{fn}: n = {n}, R = {f0:.6f}, thickness = {thick:.4f}, |grad| = {np.linalg.norm(grad):.1e}")
    print(f"   smallest Hessian eigenvalue: torus {restricted(torus):+.5f}, su(K) {restricted(suk):+.5f}, "
          f"pairs {', '.join(f'{p:+.5f}' for p in pairs)}, full {full:+.5f}")
