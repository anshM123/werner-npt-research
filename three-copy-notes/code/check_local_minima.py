"""The two equality points of the conjectured bound q_3(P_V) >= (3/4) lambda_3(V) for planes V in (C^3)^{(x)3}:
   V1 = span{GHZ, GHZ_omega}            (lambda_3 = 2/3, q_3 = 1/2),
   V2 = span{(|000> - |111>)/sqrt2, |222>} (lambda_3 = 1/2, q_3 = 3/8).
Local searches started from random perturbations (scales 1e-3 .. 0.15) in the full space of planes never go below
3/4 for the ratio q_3/lambda_3.  Requires scipy.  Run: python check_local_minima.py"""
import itertools
import numpy as np
from scipy.optimize import minimize

rng = np.random.default_rng(11)
d, D = 3, 27


def q3(X):
    T = X.reshape([d] * 6)
    tot = 0.0
    for J in itertools.product([0, 1], repeat=3):
        A, m = T, 3
        for i in reversed(range(3)):
            if J[i]:
                A = np.trace(A, axis1=i, axis2=i + m)
                m -= 1
        M = A.reshape(d ** m, d ** m)
        tot += (-0.5) ** sum(J) * np.vdot(M, M).real
    return tot


def lam3(P):
    T = P.reshape([d] * 6); out = []
    for l in range(3):
        M = T
        for i in sorted([i for i in range(3) if i != l], reverse=True):
            M = np.trace(M, axis1=i, axis2=i + M.ndim // 2)
        out.append(np.sort(np.linalg.eigvalsh(M))[-3])
    return min(out)


def plane(x):
    u = x[:D] + 1j * x[D:2 * D]; v = x[2 * D:3 * D] + 1j * x[3 * D:]
    Q, _ = np.linalg.qr(np.stack([u, v], 1))
    return Q @ Q.conj().T


def ratio(x):
    P = plane(x)
    return q3(P) / max(lam3(P), 0.05)


def ket(i, j, k):
    e = np.zeros(D, complex); e[9 * i + 3 * j + k] = 1; return e


w = np.exp(2j * np.pi / 3)
points = {"span{GHZ, GHZ_omega}": ((ket(0, 0, 0) + ket(1, 1, 1) + ket(2, 2, 2)) / np.sqrt(3),
                                   (ket(0, 0, 0) + w * ket(1, 1, 1) + w * w * ket(2, 2, 2)) / np.sqrt(3)),
          "span{(|000>-|111>)/sqrt2, |222>}": ((ket(0, 0, 0) - ket(1, 1, 1)) / np.sqrt(2), ket(2, 2, 2))}
for name, (u, v) in points.items():
    x0 = np.concatenate([u.real, u.imag, v.real, v.imag])
    best = ratio(x0)
    print(f"{name}: q_3/lambda_3 = {best:.6f}")
    for eps in (1e-3, 1e-2, 5e-2, 0.15):
        for _ in range(4):
            r = minimize(ratio, x0 + eps * rng.standard_normal(x0.size), method='L-BFGS-B', options={'maxiter': 2000})
            best = min(best, r.fun)
    print(f"   lowest value found from perturbed starts: {best:.6f}")
