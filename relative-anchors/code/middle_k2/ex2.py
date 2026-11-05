"""Exploration: soft spaces S_h = {Psi(a,h): psi_r = a_r J^-|h>, a_h = 0, a perp w}; the rest R.
Measures: F on S_h / s_h, cross S_h-S_h' / sqrt(s_h s_h'), cross S_h-R / sqrt(s_h), lambda_min(F|R)."""
import numpy as np, sys
from ex1 import reduced, kernel_basis
np.set_printoptions(linewidth=220, precision=4, suppress=True)


def soft_basis(m):
    k, n = m.k, m.np
    w = np.array(m.w, float)
    Jm = np.array(m.Jm(), float)
    blocks = []
    for h in range(k):
        # basis of {a : a_h = 0, a.w = 0}
        M = np.zeros((2, k)); M[0, h] = 1; M[1] = w
        _, _, Vt = np.linalg.svd(M)
        Ab = Vt[2:].T  # k x (k-2)
        B = np.zeros((k * n, k - 2))
        for c in range(k - 2):
            for r in range(k):
                B[r * n:(r + 1) * n, c] = Ab[r, c] * Jm[:, h]
        blocks.append(B)
    return blocks


def analyse(w, verbose=True):
    m, Y, F = reduced(list(w))
    k, n = m.k, m.np
    s = np.array(w) ** 2; S = s.sum(); x = s / S
    K = kernel_basis(m)
    Sb = soft_basis(m)
    allS = np.hstack(Sb)
    # R = orthogonal complement of span(K, allS)
    Q, _ = np.linalg.qr(np.hstack([K, allS]))
    full = np.eye(k * n)
    Rb = np.linalg.svd(full - Q @ Q.T)[0][:, :k * n - Q.shape[1]]
    # orthonormalise each soft block (Euclidean)
    Sb = [np.linalg.qr(B)[0] for B in Sb]
    FR = Rb.T @ F @ Rb
    lamR = np.linalg.eigvalsh(FR)[0]
    diag = [np.linalg.eigvalsh(B.T @ F @ B) / x[h] for h, B in enumerate(Sb)]
    cross = max(np.linalg.norm(Sb[h].T @ F @ Sb[g], 2) / np.sqrt(x[h] * x[g]) for h in range(k) for g in range(k) if g != h)
    cR = [np.linalg.norm(Sb[h].T @ F @ Rb, 2) / np.sqrt(x[h]) for h in range(k)]
    # full spectrum off kernel
    Qk, _ = np.linalg.qr(K)
    ev = np.linalg.eigvalsh((np.eye(k * n) - Qk @ Qk.T) @ F @ (np.eye(k * n) - Qk @ Qk.T) + 100 * Qk @ Qk.T)
    if verbose:
        print(f"k={k} x={np.round(x, 4)}  min eig off ker {ev[0]:.4e}  min eig/min x {ev[0] / x.min():.3f}")
        print(f"   lambda_min(F|R) = {lamR:.4f};  F|S_h / x_h eigen ranges: {[np.round(d, 3).tolist() for d in diag]}")
        print(f"   max |S_h-S_g| / sqrt(x_h x_g) = {cross:.4f};  |S_h - R| / sqrt(x_h) = {np.round(cR, 4)}")
    return ev[0], lamR, diag, cross, cR


if __name__ == '__main__':
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    for k in (4, 5):
        for trial in range(6):
            if trial == 0:
                w = np.ones(k)
            elif trial == 1:
                w = rng.random(k) + 0.3
            elif trial == 2:
                w = np.exp(rng.normal(size=k))
            elif trial == 3:
                w = np.array([1.0] * (k - 1) + [0.03])
            elif trial == 4:
                w = np.array([1.0] + [0.03] * (k - 1))
            else:
                w = np.array([1.0, 1.0] + [0.03] * (k - 2))
            analyse(w)
