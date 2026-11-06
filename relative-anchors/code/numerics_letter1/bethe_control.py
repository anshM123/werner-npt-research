"""Control for bethe_fit.py (floating point): after a small random symmetric perturbation of E_c the eigenvectors are
no longer Bethe states (the best fit leaves a residual of the size of the perturbation)."""
import numpy as np
from scipy.optimize import least_squares
from l1float import L1Block
rng = np.random.default_rng(5)
k = 6
s = np.exp(rng.normal(0, 0.7, size=k))
B = L1Block(list(s), 2)
E = B.E(2)
P = rng.normal(size=E.shape); P = (P + P.T) / 2
H2 = B.H[2]
w = np.sqrt(s)
def bethe(x):
    x1, x2 = x[0] + 1j * x[1], x[2] + 1j * x[3]
    f1 = w / (s - x1); f2 = w / (s - x2)
    return np.array([f1[a] * f2[b] + f1[b] * f2[a] for (a, b) in [tuple(sorted(y)) for y in H2]])
def fit(v):
    best = 1.0
    for trial in range(40):
        x0 = rng.normal(0, 3, 4)
        def r(x):
            b = bethe(x)
            c = np.vdot(b, v) / np.vdot(b, b)
            d = v - c * b
            return np.concatenate([d.real, d.imag])
        sol = least_squares(r, x0, max_nfev=400)
        best = min(best, np.linalg.norm(sol.fun))
    return best
for eps in [1e-2, 1e-1]:
    ev, V = np.linalg.eigh(E + eps * P)
    print(f'perturbed by {eps}: residuals', np.round([fit(V[:, n]) for n in range(0, 15, 3)], 5))
