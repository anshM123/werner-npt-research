"""Theorem 1 and Corollary 2: for a plane L containing a product (or pair) anchor, the smallest eigenvalue of q_n on
{X : ran X in L} (co-range arbitrary), relative to ||X||^2. Generic planes give a positive minimum; planes inside a
product of local qubits (a zero stratum) give 0."""
import itertools
import numpy as np
from endpoint import d, Phi

rng = np.random.default_rng(5)


def crand(*s):
    return rng.normal(size=s) + 1j * rng.normal(size=s)


def plane_min(p, u, n):
    D = d ** n
    Qb = np.linalg.qr(np.stack([p, u], 1))[0]
    basis = [np.outer(Qb[:, k], np.eye(D)[j]) for k in range(2) for j in range(D)]   # orthonormal in HS
    images = [Phi(E, n) for E in basis]
    M = np.array([[np.vdot(a, b) for b in images] for a in basis])
    return np.linalg.eigvalsh((M + M.conj().T) / 2).min()


for n in (2, 3):
    p = np.zeros(d ** n); p[0] = 1
    vals = [plane_min(p, crand(d ** n), n) for _ in range(3)]
    print(f"n={n}, product anchor |0...0>, random u: smallest eigenvalue {min(vals):.4e}")
    uS = crand(2 ** n); uS[0] = 0
    u = np.zeros(d ** n, complex)
    for idx in itertools.product(range(2), repeat=n):
        u[np.ravel_multi_index(idx, [d] * n)] = uS[np.ravel_multi_index(idx, [2] * n)]
    print(f"n={n}, product anchor, u in the local qubits span(0,1): smallest eigenvalue {plane_min(p, u, n):.1e}")
omega = np.eye(3).reshape(9) / np.sqrt(3)
for n, p in ((3, np.kron(omega, np.eye(3)[0])), (4, np.kron(omega, omega))):
    vals = [plane_min(p, crand(d ** n), n) for _ in range(2)]
    print(f"n={n}, anchor with maximally entangled pairs, random u: smallest eigenvalue {min(vals):.4e}")
