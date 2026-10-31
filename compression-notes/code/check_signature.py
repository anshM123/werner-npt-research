"""Section 7: the filter form G -> q(G_s C) has at most one negative eigenvalue.

   For an operator C and a copy s, the Hermitian form Q_s(G) = q(G_s C) on M_3 is computed in the basis of matrix
   units (G = E_ik acting on copy s from the left). The script prints its eigenvalues
   (a) for random rank-two operators at n = 2 and 3 (the endpoint holds there, so none is negative);
   (b) for the activating operator C0 of Section 5, with the form of sigma_c (x) rho_W at c = (1, 2, 2, 6.75), where
       Q(C0) < 0. At the Werner copy the form has exactly one negative eigenvalue, as Proposition 7 predicts (the pair
       state sigma_c is one-copy 2-undistillable there by the local searches of check_activation.py)."""
import itertools
import numpy as np
from werner import d, q

rng = np.random.default_rng(9)
units = [np.outer(np.eye(d)[i], np.eye(d)[k]) for i in range(d) for k in range(d)]


def filt(G, s, n):
    E = np.eye(1)
    for r in range(n):
        E = np.kron(E, G if r == s else np.eye(d))
    return E


def form_matrix(qf, C, s, n):
    """K[a, b] = polarised qf(E_a C, E_b C), with E_a the matrix units acting on copy s"""
    X = [filt(U, s, n) @ C for U in units]
    pol = lambda A, B: (qf(A + B) - qf(A - B) + 1j * qf(A - 1j * B) - 1j * qf(A + 1j * B)) / 4
    K = np.array([[pol(X[a], X[b]) for b in range(9)] for a in range(9)])
    return (K + K.conj().T) / 2


print("(a) random rank-two operators at the endpoint (no negative eigenvalue expected):")
for n in (2, 3):
    D = d ** n
    for trial in range(2):
        C = (rng.standard_normal((D, 2)) + 1j * rng.standard_normal((D, 2))) @ \
            (rng.standard_normal((2, D)) + 1j * rng.standard_normal((2, D)))
        C = C / np.linalg.norm(C)
        lows = [np.linalg.eigvalsh(form_matrix(lambda X: q(X, n, 0.5), C, s, n))[:2] for s in range(n)]
        print(f"  n = {n} #{trial}: two smallest eigenvalues per copy {[np.round(l, 5).tolist() for l in lows]}")


def gam2(c):
    return {J: sum(cS * (-1) ** sum(a * b for a, b in zip(J, S)) for S, cS in c.items()) / 4
            for J in itertools.product([0, 1], repeat=2)}


def form(X, k, coef):
    T = X.reshape([d] * (2 * k)); val = 0.0
    for J, cJ in coef.items():
        A, mm = T, k
        for i in reversed(range(k)):
            if J[i]:
                A = np.trace(A, axis1=i, axis2=i + mm); mm -= 1
        val += cJ * np.vdot(A, A).real
    return val


c = (1, 2, 2, 6.75)
coef = {J + (j3,): gJ * (-0.5) ** j3 for J, gJ in gam2(dict(zip([(0, 0), (0, 1), (1, 0), (1, 1)], c))).items()
        for j3 in (0, 1)}
Phi = np.eye(d) / np.sqrt(d)
C0 = np.einsum('ac,bd,ef->aecbfd', Phi, Phi.conj(), np.diag([1.0, 1.0, 0.0])).reshape(27, 27)
print(f"\n(b) activation example, c = {c}: Q(C0)/||C0||^2 = {form(C0, 3, coef) / np.vdot(C0, C0).real:+.6f}")
for s, name in ((2, "Werner copy"), (0, "pair 1"), (1, "pair 2")):
    ev = np.linalg.eigvalsh(form_matrix(lambda X: form(X, 3, coef), C0, s, 3))
    print(f"  {name} (copy {s}): eigenvalues {np.round(ev, 5).tolist()}  -> {int((ev < -1e-12).sum())} negative")
