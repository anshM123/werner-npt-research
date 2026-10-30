"""Crossed Cauchy-Schwarz (CCS) is not stable under tensor products.
For a Hermitian swap-invariant Z on H (x) H with F Z >= 0, say CCS(Z) holds if Tr[Z (X (x) X)] >= 0 for every Hermitian
X of rank <= 2. With Z_W = F - 1/2 (the Werner endpoint), CCS(Z_W^{(x)n}) is the n-copy endpoint for every operator whose
range and co-range lie in one plane. This script checks numerically:
  (1) the qutrit operator Z1 in ../data satisfies CCS with a positive margin, and F Z1 >= 0;
  (2) Z1 (x) Z_W (Werner factor on C^3, and also on C^2) violates CCS;
  (3) the local-unitary-invariant example Z = Z_W (x) Z_W + (3/4) I (x) I (two qutrits) satisfies CCS (it equals
      q_2(X) + (3/4)(Tr X)^2) and F Z >= 0, but Z (x) Z_W violates CCS.
Requires scipy. Run: python check_tensor_stability.py"""
import json, os
import numpy as np
from scipy.optimize import minimize

rng = np.random.default_rng(3)


def swap(d):
    F = np.zeros((d * d, d * d))
    for i in range(d):
        for j in range(d):
            F[i * d + j, j * d + i] = 1
    return F


def tensor(Z1, d1, Z2, d2):
    D = d1 * d2
    return np.kron(Z1, Z2).reshape(d1, d1, d2, d2, d1, d1, d2, d2).transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape(D * D, D * D)


def ccs_min(Z, D, starts):
    def herm2(x, s):
        u = x[:D] + 1j * x[D:2 * D]; v = x[2 * D:3 * D] + 1j * x[3 * D:]
        Y = np.outer(u, u.conj()) + s * np.outer(v, v.conj()); return Y / np.linalg.norm(Y)
    best = np.inf
    for k in range(starts):
        s = 1 if k % 2 == 0 else -1
        f = lambda x: np.trace(Z @ np.kron(herm2(x, s), herm2(x, s))).real
        best = min(best, minimize(f, rng.standard_normal(4 * D), method='L-BFGS-B', options={'maxiter': 3000}).fun)
    return best


here = os.path.dirname(os.path.abspath(__file__))
dat = json.load(open(os.path.join(here, '..', 'data', 'tensor_stability_counterexample_Z1.json')))
Z1 = np.array(dat['real']) + 1j * np.array(dat['imag'])
F3 = swap(3)
print(f"min eigenvalue of F Z1: {np.linalg.eigvalsh((F3 @ Z1 + (F3 @ Z1).conj().T) / 2).min():+.3e}")
print(f"CCS(Z1): min Tr[Z1 (X(x)X)]/||X||^2 over rank-2 Hermitian X (200 searches): {ccs_min(Z1, 3, 200):+.3e}")
for d2 in (3, 2):
    ZW = swap(d2) - 0.5 * np.eye(d2 * d2)
    print(f"CCS(Z1 (x) Z_W), Werner on C^{d2}: min = {ccs_min(tensor(Z1, 3, ZW, d2), 3 * d2, 30):+.3e}")
ZW3 = swap(3) - 0.5 * np.eye(9)
Z = tensor(ZW3, 3, ZW3, 3) + 0.75 * np.eye(81)
print(f"LU example: min eigenvalue of F Z = {np.linalg.eigvalsh(swap(9) @ Z).min():+.3e};  CCS(Z) min = {ccs_min(Z, 9, 20):+.3e}")


# For Z (x) Z_W on three qutrits, Tr[(Z (x) Z_W)(X (x) X)] = q_3(X) + (3/4) q_1(Tr_12 X) with q_1(Y) = ||Y||^2 - (Tr Y)^2/2;
# evaluate it through partial traces instead of forming the 729 x 729 operator.
def q3(X):
    T = X.reshape([3] * 6); tot = 0.0
    for J in np.ndindex(2, 2, 2):
        A, m = T, 3
        for i in reversed(range(3)):
            if J[i]:
                A = np.trace(A, axis1=i, axis2=i + m); m -= 1
        M = A.reshape(3 ** m, 3 ** m); tot += (-0.5) ** sum(J) * np.vdot(M, M).real
    return tot


def form(X):
    Y = np.trace(np.trace(X.reshape([3] * 6), axis1=0, axis2=3), axis1=0, axis2=2)
    return q3(X) + 0.75 * (np.vdot(Y, Y).real - 0.5 * np.trace(X).real ** 2)


best = np.inf
for k in range(24):
    s = 1 if k % 2 == 0 else -1
    def f(x):
        u = x[:27] + 1j * x[27:54]; v = x[54:81] + 1j * x[81:]
        X = np.outer(u, u.conj()) + s * np.outer(v, v.conj()); return form(X / np.linalg.norm(X))
    best = min(best, minimize(f, rng.standard_normal(108), method='L-BFGS-B', options={'maxiter': 3000}).fun)
print(f"LU example: CCS(Z (x) Z_W) min = {best:+.3e}")
