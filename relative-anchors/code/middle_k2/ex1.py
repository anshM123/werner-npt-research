"""Exploration: reduced form F (phi eliminated), spectrum off the kernel, soft modes vs the J^- part."""
import numpy as np, itertools, sys
from mid import Mid
np.set_printoptions(linewidth=220, precision=4, suppress=True)


def reduced(w):
    m = Mid(w, float)
    Y = np.array(m.full(), float)
    k, n = m.k, m.np
    N1 = k * n
    A = Y[:N1, :N1]; Bc = Y[:N1, N1:]; C = Y[N1:, N1:]
    # eliminate phi on w-perp (C = 2 Pi), coupling annihilates w
    F = A - Bc @ np.linalg.pinv(C) @ Bc.T
    return m, Y, F


def kernel_basis(m):
    k, n = m.k, m.np
    K = np.zeros((k * n, n))
    for i in range(n):
        for r in range(k):
            K[r * n + i, i] = m.w[r]
    return K


def jminus_part(m):
    """basis of {psi_r = J^- chi_r}: columns indexed (r, h)"""
    k, n = m.k, m.np
    Jm = np.array(m.Jm(), float)
    B = np.zeros((k * n, k * k))
    for r in range(k):
        for h in range(k):
            B[r * n:(r + 1) * n, r * k + h] = Jm[:, h]
    return B


if __name__ == '__main__':
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    for k in (4, 5):
        for trial in range(4):
            if trial == 0:
                w = np.ones(k)
            elif trial == 1:
                w = rng.random(k) + 0.3
            elif trial == 2:
                w = np.exp(rng.normal(size=k))
            else:
                w = np.array([1.0] * (k - 1) + [0.1])
            m, Y, F = reduced(list(w))
            K = kernel_basis(m)
            Qk, _ = np.linalg.qr(K)
            Pc = np.eye(len(F)) - Qk @ Qk.T
            ev, V = np.linalg.eigh(Pc @ F @ Pc + 100 * Qk @ Qk.T)
            s = w ** 2; S = s.sum()
            JB = jminus_part(m)
            Qj, _ = np.linalg.qr(JB)
            fr = [np.linalg.norm(Qj.T @ V[:, i]) ** 2 for i in range(6)]
            print(f"k={k} s/|s|={np.round(s / S, 4)}")
            print(f"   lowest eig off kernel: {ev[:8]}")
            print(f"   fraction in J^- part (lowest 6): {np.round(fr, 3)}")
            # Rayleigh on the J^- part modulo kernel
            FJ = JB.T @ F @ JB; GJ = JB.T @ JB
            # restrict to complement of kernel inside: kernel there = X = w c^T
            evj = np.sort(np.real(np.linalg.eigvals(np.linalg.pinv(GJ) @ FJ)))
            print(f"   J^- part generalized eig: {evj[:k + 4]}")
