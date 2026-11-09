"""Spectral growth (SG) and clock encodings.
For triples A, B, C, D of b x b matrices, x_{A,N}(i_1..i_N) = Tr(A_{i_1} ... A_{i_N}) and
2^N q_N(x_A x_B^T, x_C x_D^T) = Tr T(A,B;C,D)^N, T = 2 sum_ij A_i(x)B_j(x)C_i(x)D_j - sum_ij A_i(x)B_i(x)C_j(x)D_j.
A pair of rank-one operators at n copies can be written with a period-n clock (open matrix product with boundary
bond 1, made translation invariant by a cyclic shift on the bond).  The aligned clock sector of T_ab then has
cycle product 2^n h_ab, h_ab = q_n(C_a, C_b), but the sectors with misaligned clocks can have larger spectral radius.
So SG for all bond dimensions does not reduce to the n-copy Gram inequality |h_01|^2 <= h_00 h_11 this way.
usage: python sg_clock_sectors.py"""
import numpy as np
rng = np.random.default_rng(2611)


def clock_mats(x):
    """n = 2: x(i1, i2) = Tr(A_{i1} A_{i2}) summed over the two clock phases; bond C^1 + C^3."""
    X = x.reshape(3, 3); A = []
    for i in range(3):
        M = np.zeros((4, 4)); M[0, 1 + i] = 1; M[1:, 0] = X[:, i]; A.append(M)
    return A


def kron4(a, b, c, d):
    return np.kron(np.kron(np.kron(a, b), c), d)


def transfer(A, B, C, D):
    return (2 * sum(kron4(A[i], B[j], C[i], D[j]) for i in range(3) for j in range(3))
            - sum(kron4(A[i], B[i], C[j], D[j]) for i in range(3) for j in range(3)))


def Rn(C, n):
    T = C.reshape([3] * (2 * n)).astype(float)
    for s in range(n):
        Tp = np.moveaxis(T, [s, n + s], [2 * n - 2, 2 * n - 1]); tr = np.trace(Tp, axis1=-2, axis2=-1)
        T = np.moveaxis(Tp - 0.5 * tr[..., None, None] * np.eye(3), [2 * n - 2, 2 * n - 1], [s, n + s])
    return T.reshape(3 ** n, 3 ** n)


def sym_vec(x):
    return (x.reshape(3, 3) + x.reshape(3, 3).T).reshape(9)


worst = 0.0
aligned_ok = True
for trial in range(200):
    xs = [rng.normal(size=9) for _ in range(4)]
    M = [clock_mats(x) for x in xs]
    # identity check: the clock MPS produces the rotation-symmetrised vectors
    ok = all(np.allclose([np.trace(M[k][i] @ M[k][j]) for i in range(3) for j in range(3)], sym_vec(xs[k])) for k in range(4))
    T00 = transfer(M[0], M[1], M[0], M[1]); T11 = transfer(M[2], M[3], M[2], M[3])
    C0, C1 = np.outer(xs[0], xs[1]), np.outer(xs[2], xs[3])
    h00, h11 = np.sum(C0 * Rn(C0, 2)), np.sum(C1 * Rn(C1, 2))
    r00 = np.abs(np.linalg.eigvals(T00)).max() ** 2; r11 = np.abs(np.linalg.eigvals(T11)).max() ** 2
    # the aligned clock sector contributes the eigenvalue 4 h_aa of T_aa^2; misaligned sectors may exceed it
    al = min(np.abs(np.linalg.eigvals(T00 @ T00) - 4 * h00).min() / (4 * h00), np.abs(np.linalg.eigvals(T11 @ T11) - 4 * h11).min() / (4 * h11))
    aligned_ok = aligned_ok and al < 1e-9
    worst = max(worst, r00 / (4 * h00), r11 / (4 * h11))
    if trial < 3:
        print(f"trial {trial}: clock vectors = symmetrised vectors: {ok}; rho(T_00)^2 = {r00:.3f} vs 4 h_00 = {4*h00:.3f};"
              f" rho(T_11)^2 = {r11:.3f} vs 4 h_11 = {4*h11:.3f}")
print(f"4 h_aa is an eigenvalue of T_aa^2 in every trial: {aligned_ok}")
print(f"largest ratio rho(T_aa)^2 / (4 h_aa) over 200 random two-copy pairs: {worst:.3f} (> 1: a misaligned sector dominates)")
