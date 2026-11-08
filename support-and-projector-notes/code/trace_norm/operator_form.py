"""Section 9, the operator form (9.8).  For unit psi in C^2 (x) (C^3)^n, with Phi'(Z) = Tr(Z) I - tZ, R(Z) = Tr(Z) I - Z,
tau = t/(1-t) and rho_RK the marginal of psi psi^* on the reference qubit and the copies in K, check the identity
  t (1-t)^{n-1} [ sum_{K != 0} tau^{|K|-1} (id (x) R^{(x)K})(rho_RK) (x) I + (I_2 - rho_R) (x) I ]
      = (id (x) Phi'^{(x)n})(psi psi^*) - (1-t)^{n-1} (rho_R - t I_2) (x) I
for n = 1, 2, 3 and t = 0.3, 0.4, 0.5 on random psi (floating point).
usage: python operator_form.py"""
import itertools, numpy as np
rng = np.random.default_rng(2)
def omega_direct(psi, n, t):
    m = 1 - t; D = 3 ** n
    X = np.outer(psi, psi.conj()).reshape([2] + [3] * n + [2] + [3] * n)
    # apply Phi' = Tr(.)I - t id at each site
    Y = X.copy()
    for s in range(n):
        ax_r, ax_c = 1 + s, 2 + n + s
        tr = np.trace(Y, axis1=ax_r, axis2=ax_c)
        Y = -t * Y + np.moveaxis(np.multiply.outer(tr, np.eye(3)), [-2, -1], [ax_r, ax_c])
    L = Y.reshape(2 * D, 2 * D)
    rhoR = np.einsum('aibi->ab', np.outer(psi, psi.conj()).reshape(2, D, 2, D))
    return L - m ** (n - 1) * np.kron(rhoR - t * np.eye(2), np.eye(D))
def omega_R(psi, n, t):
    m = 1 - t; tau = t / m; D = 3 ** n
    rho = np.outer(psi, psi.conj()).reshape([2] + [3] * n + [2] + [3] * n)
    tot = np.zeros((2 * D, 2 * D), complex)
    for r in range(1, n + 1):
        for K in itertools.combinations(range(n), r):
            Y = rho.copy()
            for s in range(n):
                ax_r, ax_c = 1 + s, 2 + n + s
                tr = np.trace(Y, axis1=ax_r, axis2=ax_c)
                full = np.moveaxis(np.multiply.outer(tr, np.eye(3)), [-2, -1], [ax_r, ax_c])
                Y = (full - Y) if s in K else full          # R = Tr(.)I - id on K, Tr(.)I elsewhere (gives rho_RK (x) I)
            tot += tau ** (r - 1) * Y.reshape(2 * D, 2 * D)
    rhoR = np.einsum('aibi->ab', np.outer(psi, psi.conj()).reshape(2, D, 2, D))
    return tot + np.kron(np.eye(2) - rhoR, np.eye(D))
worst = 0
for n in (1, 2, 3):
    for t in (0.3, 0.4, 0.5):
        m = 1 - t
        for trial in range(3):
            psi = rng.normal(size=2 * 3 ** n) + 1j * rng.normal(size=2 * 3 ** n); psi /= np.linalg.norm(psi)
            A = omega_direct(psi, n, t); B = omega_R(psi, n, t)
            worst = max(worst, np.abs(A - t * m ** (n - 1) * B).max())
print("identity behind (9.8), n = 1, 2, 3: max error", worst)
print("PASS" if worst < 1e-12 else "FAIL")
