"""Section 3: the projector case as a monogamy relation, every n.

   For a plane V with orthonormal basis v0, v1 put psi = (|0>_R v0 + |1>_R v1)/sqrt2, with R a qubit, and
   c_J = || psi_{JR} - psi_J (x) I/2 ||_2^2 for subsets J of the copies. Then
       q_n(1/2, P_V) = (16/3) sum_J [ (-1/2)^{|J|} + (1/2) (-1/2)^{n-|J|} ] c_J.
   Checked for random planes, n = 1..4 (d = 3) and n = 1..3 (d = 4)."""
import itertools
import numpy as np

rng = np.random.default_rng(5)


def q_direct(X, n, dd, t=0.5):
    T = X.reshape([dd] * (2 * n)); tot = 0.0
    for J in itertools.product([0, 1], repeat=n):
        A, mm = T, n
        for i in reversed(range(n)):
            if J[i]:
                A = np.trace(A, axis1=i, axis2=i + mm); mm -= 1
        tot += (-t) ** sum(J) * np.vdot(A, A).real
    return tot


def monogamy_form(V, n, dd):
    psi = np.stack([V[:, 0], V[:, 1]]) / np.sqrt(2)        # index order: R, copies
    T = psi.reshape([2] + [dd] * n)
    total = 0.0
    for J in itertools.product([0, 1], repeat=n):
        traced = [1 + i for i in range(n) if not J[i]]
        rho = np.tensordot(T, T.conj(), axes=(traced, traced))
        dJ = dd ** sum(J)
        rho = rho.reshape(2 * dJ, 2 * dJ)
        rhoJ = np.einsum('aiaj->ij', rho.reshape(2, dJ, 2, dJ))
        c = np.linalg.norm(rho - np.kron(np.eye(2) / 2, rhoJ)) ** 2
        total += ((-0.5) ** sum(J) + 0.5 * (-0.5) ** (n - sum(J))) * c
    return 16 / 3 * total


for dd, nmax in ((3, 4), (4, 3)):
    for n in range(1, nmax + 1):
        D = dd ** n
        V = np.linalg.qr(rng.standard_normal((D, 2)) + 1j * rng.standard_normal((D, 2)))[0]
        print(f"d={dd} n={n}: q_n(P_V) = {q_direct(V @ V.conj().T, n, dd):.10f}   monogamy form = {monogamy_form(V, n, dd):.10f}")
