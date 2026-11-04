"""Shared helpers: the endpoint form on single-excitation and tangent ranges, written through qubit POVMs.

Range vectors u_r = a_r e + sum_s (U_s e_r)@s, e = |0...0>, U_s a 2x2 block into E_s = span{|1>, |2>} at copy s.
In the gauge |abar><abar| + sum_s U_s^* U_s = I the reference Gram matrix is the identity, and
    m = min over co-ranges of q_n(C)/||C||^2 = lambda_min K,   K = [Psi^(x)n(u_r u_r'^*)]_{r, r'},
with Psi(Y) = Tr(Y) I - Y/2 on each copy (the 2-positivity form).  Up to local unitaries on the E_s, m depends only on
the POVM (A_0 = |abar><abar|, A_s = U_s^* U_s)."""
import numpy as np, itertools

def strings(n):
    return list(itertools.product(range(3), repeat=n))

def psi_n(X, n):
    """Psi^{(x)n}(X), Psi(Y) = Tr(Y) I - Y/2 on each copy."""
    T = X.reshape([3] * (2 * n))
    for s in range(n):
        tr = np.trace(T, axis1=s, axis2=n + s)
        full = np.moveaxis(np.einsum('ab,...->...ab', np.eye(3), tr), [2 * n - 2, 2 * n - 1], [s, n + s])
        T = full - 0.5 * T
    return T.reshape(3 ** n, 3 ** n)

def lam_map(X, n, t=0.5):
    """Lambda^{(x)n}(X), Lambda(Y) = Y - t Tr(Y) I; q_n(X) = <X, Lambda^{(x)n}(X)>."""
    T = X.reshape([3] * (2 * n))
    for s in range(n):
        tr = np.trace(T, axis1=s, axis2=n + s)
        full = np.moveaxis(np.einsum('ab,...->...ab', np.eye(3), tr), [2 * n - 2, 2 * n - 1], [s, n + s])
        T = T - t * full
    return T.reshape(3 ** n, 3 ** n)

def q(X, n):
    return np.vdot(X, lam_map(X, n)).real

def sqrt_psd(a):
    w, V = np.linalg.eigh(a); return V @ np.diag(np.sqrt(np.clip(w, 0, None))) @ V.conj().T

def tangent_range(U, a, n):
    idx = {s: i for i, s in enumerate(strings(n))}; u = np.zeros((2, 3 ** n), complex)
    for r in range(2):
        u[r, 0] += a[r]
        for s in range(n):
            for z in (1, 2):
                t = [0] * n; t[s] = z; u[r, idx[tuple(t)]] += U[s][z - 1, r]
    return u

def K_of(u, n):
    D = 3 ** n; K = np.zeros((2 * D, 2 * D), complex)
    for r in range(2):
        for rp in range(2):
            K[r * D:(r + 1) * D, rp * D:(rp + 1) * D] = psi_n(np.outer(u[r], u[rp].conj()), n)
    return (K + K.conj().T) / 2

def m_povm(A, n, a=None):
    """m for the POVM (A_1..A_n) (single excitations, a = None) or with vacuum amplitude a (tangent class);
    the gauge must hold: sum A_s (+ abar abar^*) = I."""
    U = [sqrt_psd(x) for x in A]; u = tangent_range(U, np.zeros(2) if a is None else a, n)
    return np.linalg.eigvalsh(K_of(u, n))[0]

def sector_mu(A, n, c):
    """lambda_min of K restricted to co-range strings with c excitations (single-excitation range)."""
    D = 3 ** n; u = tangent_range([sqrt_psd(x) for x in A], np.zeros(2), n); K = K_of(u, n)
    idx = np.where(np.array([sum(x != 0 for x in s) for s in strings(n)]) == c)[0]; sel = np.concatenate([idx, D + idx])
    return np.linalg.eigvalsh(K[np.ix_(sel, sel)])[0]

def povm_from(x, n):
    """Interior qubit POVM from unconstrained real parameters x (length 8n)."""
    G = (x[:4 * n] + 1j * x[4 * n:]).reshape(n, 2, 2); A = [g @ g.conj().T for g in G]; S = sum(A)
    w, V = np.linalg.eigh(S); W = V @ np.diag(w ** -0.5) @ V.conj().T; return [W @ x_ @ W for x_ in A]

def tangent_from(x, n):
    """(U, a) in the gauge from unconstrained parameters x (length 2(4n+2))."""
    z = x[:len(x) // 2] + 1j * x[len(x) // 2:]; U = [z[4 * s:4 * s + 4].reshape(2, 2) for s in range(n)]; a = z[4 * n:4 * n + 2]
    S = np.outer(a.conj(), a) + sum(u.conj().T @ u for u in U); w, V = np.linalg.eigh(S); T = V @ np.diag(w ** -0.5) @ V.conj().T
    return [u @ T for u in U], T.T @ a

def m_tangent(x, n):
    U, a = tangent_from(x, n); return np.linalg.eigvalsh(K_of(tangent_range(U, a, n), n))[0]
