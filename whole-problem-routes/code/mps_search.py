"""Check for Observation 7.4 (whole-problem-routes/README.md): translation-invariant real MPS inputs x, y (bond dimension D, periodic) for the map
P_a(X) = a Tr(X) I - (X + X^T)/4 on n qutrit copies.  Value R = <y|P^n(xx^T)|y> / (|x|^2 |y|^2) = Tr(T^n)/(Tr(Ex^n) Tr(Ey^n)),
T = sum_{ijkl} P[k,l;i,j] A_i (x) A_j (x) B_k (x) B_l,  P[k,l;i,j] = <e_k|P(e_i e_j^T)|e_l> = a d_ij d_kl - (d_ki d_lj + d_kj d_li)/4.
Minimises R by BFGS from random starts and from the plane witness (GHZ in the Y basis: A_0 = I/sqrt2,
A_1 = [[0,1],[-1,0]]/sqrt2, A_2 = 0, padded to bond dimension D).
usage: python mps_search.py n D starts a1 a2 ...   (logs/mps_<n>_<D>.log; arguments in the first line of each log)"""
import sys
import numpy as np
from scipy.optimize import minimize

n, D, starts = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
avals = [float(s) for s in sys.argv[4:]]
rng = np.random.default_rng(1000 * n + D)


def Pten(a):
    P = np.zeros((3, 3, 3, 3))   # P[k, l, i, j]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    P[k, l, i, j] = a * (i == j) * (k == l) - ((k == i) * (l == j) + (k == j) * (l == i)) / 4
    return P


def value(params, P):
    A = params[:3 * D * D].reshape(3, D, D); B = params[3 * D * D:].reshape(3, D, D)
    Ex = sum(np.kron(A[i], A[i]) for i in range(3)); Ey = sum(np.kron(B[i], B[i]) for i in range(3))
    sx = max(np.max(np.abs(np.linalg.eigvals(Ex))), 1e-300); sy = max(np.max(np.abs(np.linalg.eigvals(Ey))), 1e-300)
    A = A / np.sqrt(sx); B = B / np.sqrt(sy)
    Ex = Ex / sx; Ey = Ey / sy
    T = np.zeros((D ** 4, D ** 4))
    for i in range(3):
        for j in range(3):
            AA = np.kron(A[i], A[j])
            for k in range(3):
                for l in range(3):
                    if P[k, l, i, j] != 0:
                        T += P[k, l, i, j] * np.kron(AA, np.kron(B[k], B[l]))
    nx = np.trace(np.linalg.matrix_power(Ex, n)); ny = np.trace(np.linalg.matrix_power(Ey, n))
    if nx < 1e-12 or ny < 1e-12:
        return 1.0
    return np.trace(np.linalg.matrix_power(T, n)) / (nx * ny)


def ghz_seed():
    A = np.zeros((3, D, D))
    A[0, :2, :2] = np.eye(2) / np.sqrt(2); A[1, :2, :2] = np.array([[0, 1], [-1, 0]]) / np.sqrt(2)
    return A


for a in avals:
    P = Pten(a)
    best = np.inf; seeds = []
    g = ghz_seed()
    seeds.append(np.concatenate([g.ravel(), g.ravel()]))
    for _ in range(starts):
        seeds.append(rng.normal(size=2 * 3 * D * D))
    for s0 in seeds:
        r = minimize(value, s0, args=(P,), method='BFGS', options={'maxiter': 400, 'gtol': 1e-10})
        best = min(best, r.fun)
    an = (1 + 2 ** (1 - 1 / n)) / 4
    ghz = 2.0 ** -n * ((2 * a - 0.5) ** n - 0.5)
    print(f'n = {n}, D = {D}, a = {a:.5f} (a_n = {an:.5f}): least value found {best:+.4e}; plane witness value {ghz:+.4e}; '
          f'product value (a - 1/2)^n = {(a - 0.5) ** n:.4e}', flush=True)
