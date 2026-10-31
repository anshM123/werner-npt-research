"""Section 6: sums of two block-product operators (a Schur product argument).

   Partition the copies into blocks of size <= 2 and take C = sum_{a=1,2} (x)_j |alpha_a^(j)><beta_a^(j)| with
   arbitrary (entangled) vectors on each block. Then q_n(1/2, C) = 1^T (Hadamard product of the 2x2 Gram matrices
   N^(j)) 1, each N^(j) is positive semidefinite by the two-copy theorem, so q_n >= 0. The script checks the
   identity, prints q_n/||C||^2, the smallest eigenvalue of each N^(j), and the local spectra (thickness) at every
   copy, for n = 4 (blocks 2+2) and n = 5 (blocks 2+2+1)."""
import numpy as np
from werner import d, q, q_pair

rng = np.random.default_rng(8)
t = 0.5


def rvec(k):
    v = rng.standard_normal(d ** k) + 1j * rng.standard_normal(d ** k)
    return v / np.linalg.norm(v)


def local_spectra(C, n):
    out = []
    M = C @ C.conj().T + C.conj().T @ C
    for s in range(n):
        T = M.reshape([d] * (2 * n)); T = np.moveaxis(T, [s, n + s], [0, n]).reshape(d, d ** (n - 1), d, d ** (n - 1))
        ev = np.linalg.eigvalsh(np.einsum('aibi->ab', T))
        out.append(round(ev[0] / ev.sum(), 4))
    return out


for blocks in ([2, 2], [2, 2, 1]):
    n = sum(blocks)
    for trial in range(3):
        facs = [[(rvec(k), rvec(k)) for k in blocks] for a in range(2)]          # facs[a][j] = (alpha, beta)
        terms = []
        for a in range(2):
            T = np.eye(1)
            for (al, be) in facs[a]:
                T = np.kron(T, np.outer(al, be.conj()))
            terms.append(T)
        C = terms[0] + terms[1]
        Ns = []
        for j, k in enumerate(blocks):
            A = [np.outer(facs[a][j][0], facs[a][j][1].conj()) for a in range(2)]
            Ns.append(np.array([[q_pair(A[a], A[b], k, t) for b in range(2)] for a in range(2)]))
        had = Ns[0]
        for N in Ns[1:]:
            had = had * N
        print(f"blocks {blocks} #{trial}: q_n = {q(C, n, t):+.8f}, Schur form = {had.sum().real:+.8f}, "
              f"q_n/||C||^2 = {q(C, n, t) / np.vdot(C, C).real:+.4f}, min eig of the N^(j): "
              f"{[round(np.linalg.eigvalsh((N + N.conj().T) / 2)[0], 6) for N in Ns]}, "
              f"smallest local eigenvalue / trace per copy: {local_spectra(C, n)}", flush=True)

# general form: range of C inside (x)_j V_j (V_j a plane in each two-copy block), arbitrary co-range, any rank
print("\nrange in a product of block planes, arbitrary co-range (rank up to 2^m):")
for blocks in ([2, 2], [2, 2, 1]):
    n = sum(blocks)
    for trial in range(3):
        P = np.eye(1)
        for k in blocks:
            V = np.linalg.qr(rng.standard_normal((d ** k, 2)) + 1j * rng.standard_normal((d ** k, 2)))[0]
            P = np.kron(P, V @ V.conj().T)
        G = rng.standard_normal((d ** n, d ** n)) + 1j * rng.standard_normal((d ** n, d ** n))
        C = P @ G                                     # range inside (x)_j V_j, generic co-range
        print(f"blocks {blocks} #{trial}: rank {np.linalg.matrix_rank(C)}, q_n/||C||^2 = {q(C, n, t) / np.vdot(C, C).real:+.5f}",
              flush=True)
