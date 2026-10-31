"""Section 2: the zeros C0 = P (x) |a><b| at t = 1/2 (P = |0><0| + |1><1| on the first site, a, b any unit vectors on the
   other n-1 sites).  f(Y, Z) = q_n(Y Z^dagger) with Y = [|0>a, |1>a], Z = [|0>b, |1>b] has zero gradient there; we print
   the smallest eigenvalues of its Hessian (finite differences of the exact gradient) for random a, b, n = 3, 4.
   Proposition 4 says they are >= 0 (up to rounding) whenever q_{n-1} >= 0 on rank two, which is known for n - 1 <= 3."""
import sys
import numpy as np
from werner import d, applyB, unpack, pack

t = 0.5
rng = np.random.default_rng(1)
NR = int(sys.argv[1]) if len(sys.argv) > 1 else 3


def grad(x, n):
    D = d ** n
    Y, Z = unpack(x, D); C = Y @ Z.conj().T
    g = applyB(C, n, t)
    gY = g @ Z; gZ = g.conj().T @ Y
    return 2 * np.concatenate([gY.real.ravel(), gY.imag.ravel(), gZ.real.ravel(), gZ.imag.ravel()])


for n in (3, 4):
    E = d ** (n - 1)
    for trial in range(NR):
        a = rng.standard_normal(E) + 1j * rng.standard_normal(E); a /= np.linalg.norm(a)
        b = rng.standard_normal(E) + 1j * rng.standard_normal(E); b /= np.linalg.norm(b)
        Y = np.stack([np.kron(np.eye(d)[k], a) for k in range(2)], axis=1)
        Z = np.stack([np.kron(np.eye(d)[k], b) for k in range(2)], axis=1)
        x0 = pack(Y, Z); N = x0.size; h = 1e-6
        g0 = np.linalg.norm(grad(x0, n))
        H = np.zeros((N, N))
        for i in range(N):
            e = np.zeros(N); e[i] = h
            H[:, i] = (grad(x0 + e, n) - grad(x0 - e, n)) / (2 * h)
        ev = np.linalg.eigvalsh((H + H.T) / 2)
        print(f"n={n}, random entangled a, b #{trial}: |gradient| = {g0:.1e}; smallest Hessian eigenvalues "
              f"{np.round(ev[:4], 9)}; number below -1e-8: {(ev < -1e-8).sum()}", flush=True)
