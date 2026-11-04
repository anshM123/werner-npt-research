"""Shared helpers: the endpoint form q_n, tangent-class operators, single-copy squeeze forms."""
import numpy as np, itertools
d = 3


def strings(n):
    return list(itertools.product(range(d), repeat=n))


def Lmap(X, n, t=0.5):
    """Lambda^{(x)n}(X), Lambda(Y) = Y - t Tr(Y) I on each copy; q_n(X) = <X, Lambda^{(x)n}(X)>."""
    T = X.reshape([d] * (2 * n))
    for s in range(n):
        tr = np.trace(T, axis1=s, axis2=n + s)
        full = np.moveaxis(np.einsum('ab,...->...ab', np.eye(d), tr), [2 * n - 2, 2 * n - 1], [s, n + s])
        T = T - t * full
    return T.reshape(d ** n, d ** n)


def q(X, n):
    return np.vdot(X, Lmap(X, n)).real


def q_direct(X, n):
    """sum_J (-1/2)^|J| ||Tr_J X||^2, by explicit partial traces (independent of Lmap)."""
    T = X.reshape([d] * (2 * n)); tot = 0.0
    for J in itertools.product([0, 1], repeat=n):
        Y = T; cur = list(range(n))
        for s in sorted([s for s in range(n) if J[s]], reverse=True):
            i = cur.index(s); m = len(cur); Y = np.trace(Y, axis1=i, axis2=m + i); cur.pop(i)
        tot += (-0.5) ** sum(J) * np.sum(abs(Y) ** 2)
    return tot


def tangent_range(U, a, n):
    """Range vectors u_r = a_r |0..0> + sum_s (U_s e_r) @ copy s, U_s: 2x2 (rows = excitation letters 1, 2)."""
    st = strings(n); idx = {s: i for i, s in enumerate(st)}; D = d ** n
    u = np.zeros((2, D), complex)
    for r in range(2):
        u[r, 0] += a[r]
        for s in range(n):
            for z in (1, 2):
                t = [0] * n; t[s] = z; u[r, idx[tuple(t)]] += U[s][z - 1, r]
    return u


# Herm_3 basis: diagonals, then (0,1), (0,2), (1,2) real and imaginary off-diagonals
B = []
for i in range(3):
    E = np.zeros((3, 3), complex); E[i, i] = 1; B.append(E)
for i in range(3):
    for j in range(i + 1, 3):
        E = np.zeros((3, 3), complex); E[i, j] = E[j, i] = 1 / np.sqrt(2); B.append(E)
        E = np.zeros((3, 3), complex); E[i, j] = -1j / np.sqrt(2); E[j, i] = 1j / np.sqrt(2); B.append(E)
Bt = np.array(B)


def _ptrace_keep(X, n, J):
    T = X.reshape([d] * (2 * n)); cur = list(range(n))
    for t in sorted(J, reverse=True):
        i = cur.index(t); m = len(cur); T = np.trace(T, axis1=i, axis2=m + i); cur.pop(i)
    return T, cur


def squeeze_forms(X, n, s):
    """Quadratic forms M -> q_n(G X G), M -> ||G X G||^2 (G = M^{1/2} on copy s) as 9x9 matrices in the basis Bt:
    Q(M) = sum_{J not s} (-1/2)^|J| Tr(Y_J^* M Y_J M) + sum_{J containing s} (-1/2)^|J| ||Tr_J(X M)||^2, N(M) = Tr(X^* M X M)."""
    K = np.zeros((9, 9)); Kn = np.zeros((9, 9))
    others = [t for t in range(n) if t != s]
    for Jbits in itertools.product([0, 1], repeat=n - 1):
        J = [others[i] for i in range(n - 1) if Jbits[i]]; c = (-0.5) ** len(J)
        T, cur = _ptrace_keep(X, n, J); k = len(cur); i = cur.index(s)
        T = np.moveaxis(T, [i, k + i], [0, 1]).reshape(3, 3, -1)
        G = np.einsum('abr,cdr->acbd', T.conj(), T)
        val = np.einsum('kac,lfb,acbf->kl', Bt, Bt, G).real
        K += c * val
        if not J:
            Kn += val
        Z = np.einsum('kba,abr->kr', Bt, T)
        K += (-0.5 * c) * np.real(Z.conj() @ Z.T)
    return (K + K.T) / 2, (Kn + Kn.T) / 2


def pencil_margin(K, Kn, lo=-5.0, hi=5.0):
    """nu = max_lam lambda_min(K - lam Kn) / ||Kn||; nu < 0 means K - lam Kn is PSD for no lam in [lo, hi]."""
    from scipy.optimize import minimize_scalar
    sc = np.linalg.norm(Kn, 2)
    r = minimize_scalar(lambda l: -np.linalg.eigvalsh(K - l * Kn)[0], bounds=(lo, hi), method='bounded', options={'xatol': 1e-12})
    return -r.fun / sc, r.x
