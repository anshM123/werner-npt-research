"""The segment between the Werner endpoint and the critical map (whole-problem-routes/README.md, Section 5).
On C^3 (x) C^3 let P_0, P_1, P_2 project onto the maximally entangled vector, the antisymmetric and the traceless
symmetric vectors, and X_w = -P_0 + w P_1 + 2 P_2 = (1 + w/2) I + (1 - w/2) F - Omega.
  w = 2: twice the Werner endpoint operator I - Omega/2;  w = 3: (12/2) Choi(T.F) of the 5/6 completion;
  w = 4: four times the Choi matrix of P_{3/4}.
 (a) least value of <vec C| X_w^(x)n |vec C> / ||C||^2 over rank-two C, n = 1, 2, 3 (and 4 with the flag --n4);
 (b) 2 P_{3/4} = Lambda + W with Lambda = D - id/2 and W = (D - T)/2; least output eigenvalue of Lambda^(x)a (x) W^(x)b.
usage: python segment_check.py [--n4]"""
import sys
import numpy as np
from scipy.optimize import minimize

rng = np.random.default_rng(424)
I3 = np.eye(3)


def site(T, n, s, f):
    Tp = np.moveaxis(T, [s, n + s], [2 * n - 2, 2 * n - 1])
    return np.moveaxis(f(Tp), [2 * n - 2, 2 * n - 1], [s, n + s])


def lin(a, b, c):
    """X -> a Tr(X) I + b X + c X^T on the last two axes."""
    return lambda Tp: b * Tp + c * np.swapaxes(Tp, -1, -2) + a * np.trace(Tp, axis1=-2, axis2=-1)[..., None, None] * I3


def apply(X, n, maps):
    T = X.reshape([3] * (2 * n))
    for s in range(n): T = site(T, n, s, maps[s])
    return T.reshape(3 ** n, 3 ** n)


def least_rank_two(n, w, starts):
    N = 3 ** n; M = [lin(-1.0, 1 + w / 2, 1 - w / 2)] * n; best = np.inf

    def f(x):
        U = (x[:2 * N] + 1j * x[2 * N:4 * N]).reshape(N, 2); V = (x[4 * N:6 * N] + 1j * x[6 * N:]).reshape(N, 2)
        C = U @ V.conj().T; nr = np.vdot(C, C).real; G = apply(C, n, M); val = np.vdot(C, G).real
        Gr = (G - val / nr * C) / nr; gU = 2 * Gr @ V; gV = 2 * Gr.conj().T @ U
        return val / nr, np.concatenate([gU.real.ravel(), gU.imag.ravel(), gV.real.ravel(), gV.imag.ravel()])
    for s in range(starts):
        x0 = rng.normal(size=8 * N)
        if s % 2: x0[2 * N:4 * N] = 0; x0[6 * N:] = 0
        best = min(best, minimize(f, x0, jac=True, method='L-BFGS-B', options={'maxiter': 4000, 'gtol': 1e-11, 'ftol': 1e-15}).fun)
    return best


def least_output(maps, starts):
    n = len(maps); D = 3 ** n; best = np.inf
    for s in range(starts):
        psi = rng.normal(size=D) + 1j * rng.normal(size=D)
        if s % 3 == 1: psi = psi.real + 0j
        psi /= np.linalg.norm(psi); old = np.inf
        for _ in range(400):
            Y = apply(np.outer(psi, psi.conj()), n, maps); e1, U1 = np.linalg.eigh((Y + Y.conj().T) / 2); phi = U1[:, 0]
            Z = apply(np.outer(phi, phi.conj()), n, maps); e2, U2 = np.linalg.eigh((Z + Z.conj().T) / 2); psi = U2[:, 0]
            if abs(e2[0] - old) < 1e-15: break
            old = e2[0]
        best = min(best, e2[0])
    return best


print("(a) least value of the X_w form on rank-two C")
for n, starts in ((1, 6), (2, 12), (3, 12)) + (((4, 6),) if '--n4' in sys.argv else ()):
    ws = (1.8, 2.0, 3.0, 4.0, 4.2, 5.0) if n < 4 else (2.0, 3.0, 4.0, 4.2)
    print(f"  n = {n}: " + ", ".join(f"w = {w}: {least_rank_two(n, w, starts):+.5f}" for w in ws), flush=True)
Lam = lin(1.0, -0.5, 0.0); W = lin(0.5, 0.0, -0.5); P34 = lin(0.75, -0.25, -0.25)
X = (rng.normal(size=(1, 1, 3, 3)) + 1j * rng.normal(size=(1, 1, 3, 3)))
kraus = sum((np.outer(I3[i], I3[j]) - np.outer(I3[j], I3[i])) @ X[0, 0] @ (np.outer(I3[i], I3[j]) - np.outer(I3[j], I3[i])).T
            for i in range(3) for j in range(i + 1, 3)) / 2
print(f"(b) 2 P_(3/4) = Lambda + W: {np.allclose(2 * P34(X), Lam(X) + W(X))};"
      f" W(X) = (1/2) sum_(i<j) A_ij X A_ij^T, A_ij = E_ij - E_ji: {np.allclose(kraus, W(X)[0, 0])}")
for a, b in ((1, 1), (2, 1), (3, 1), (1, 2), (2, 2), (1, 3)):
    print(f"  Lambda^{a} (x) W^{b}: least output eigenvalue {least_output([Lam] * a + [W] * b, 12 if a + b < 4 else 8):+.6f}", flush=True)
