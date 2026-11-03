"""Stage 3f: induction16 (5) balanced-factor sum-of-squares identity for the two-factor flat margin."""
import numpy as np
rng = np.random.default_rng(7)
PA = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])]


def ptr(W, d1, d2, which):
    T = W.reshape(d1, d2, d1, d2)
    return np.einsum('ajbj->ab', T) if which == 2 else np.einsum('iaib->ab', T)


def Sf(W, d1, d2):
    return np.linalg.norm(ptr(W, d1, d2, 1)) ** 2 + np.linalg.norm(ptr(W, d1, d2, 2)) ** 2


def hf(W, d1, d2):
    return np.linalg.norm(W) ** 2 - Sf(W, d1, d2) + abs(np.trace(W)) ** 2


worst_id, worst_ineq, worst_h = 0, 0, 0
for t in range(400):
    d1, d2 = rng.integers(2, 5), rng.integers(2, 5)
    n = d1 * d2
    Z = (rng.normal(size=(n, 2)) + 1j * rng.normal(size=(n, 2))) @ (rng.normal(size=(2, n)) + 1j * rng.normal(size=(2, n)))
    U, s, Vh = np.linalg.svd(Z)
    X0 = U[:, :2] * np.sqrt(s[:2]); Y0 = Vh.conj().T[:, :2] * np.sqrt(s[:2])
    M = Y0.conj().T @ X0
    m = np.array([np.trace(P @ M) / 2 for P in PA])           # M = Tr(M)/2 I + m.sigma
    nvec = np.cross(m.real, m.imag)
    if np.linalg.norm(nvec) < 1e-12:
        nvec = np.array([1.0, 0, 0]) if abs(m.real[0]) < 0.9 * np.linalg.norm(m.real) + 1e-30 else np.array([0, 1.0, 0])
        nvec = nvec - (nvec @ m.real) * m.real / max(m.real @ m.real, 1e-30)
    nvec = nvec / np.linalg.norm(nvec)
    w, Q = np.linalg.eigh(sum(nvec[k] * PA[k] for k in range(3)))
    X = X0 @ Q; Y = Y0 @ Q
    x1, x2, y1, y2 = X[:, 0], X[:, 1], Y[:, 0], Y[:, 1]
    G = X.conj().T @ X
    assert np.allclose(G, Y.conj().T @ Y) and np.isclose(np.vdot(y1, x1), np.trace(Z) / 2) and np.isclose(np.vdot(y2, x2), np.trace(Z) / 2)
    nuc = s[0] + s[1]
    lhs = np.linalg.norm(Z) ** 2 + nuc ** 2 / 2 + abs(np.trace(Z)) ** 2 / 2 - Sf(Z, d1, d2)
    o = np.outer
    rhs = (G[0, 0].real - G[1, 1].real) ** 2 / 2 + hf(o(x1, y1.conj()), d1, d2) + hf(o(x2, y2.conj()), d1, d2) \
        + hf(o(x2, x1.conj()), d1, d2) + hf(o(y2, y1.conj()), d1, d2) + Sf(o(x2, x1.conj()) - o(y2, y1.conj()), d1, d2)
    worst_id = max(worst_id, abs(lhs - rhs) / np.linalg.norm(Z) ** 2)
    q2 = np.linalg.norm(Z) ** 2 - 0.5 * Sf(Z, d1, d2) + 0.25 * abs(np.trace(Z)) ** 2
    worst_ineq = max(worst_ineq, abs((2 * q2 - (s[0] - s[1]) ** 2 / 2) - lhs) / np.linalg.norm(Z) ** 2)
    worst_h = min(worst_h, min(hf(o(a, b.conj()), d1, d2) for a, b in ((x1, y1), (x2, y2), (x2, x1), (y2, y1))) / np.linalg.norm(Z) ** 2)
msg = (f"induction16 (5): max |LHS-RHS|/||Z||^2 = {worst_id:.2e}; max |LHS - (2q2-(s1-s2)^2/2)|/||Z||^2 = {worst_ineq:.2e}; "
       f"min h(rank-one)/||Z||^2 = {worst_h:.2e} (>= 0)  [400 random rank-2 Z, d1,d2 in 2..4]")
print(msg); open('stage3f_sos.log', 'w').write(msg + '\n')
