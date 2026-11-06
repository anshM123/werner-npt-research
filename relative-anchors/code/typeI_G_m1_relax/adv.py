"""Exploration (floating point): adversarial search for small eigenvalues of the relaxed G_rel (form Psi of NOTES.md)
in the frame of Section 1.  sigma = 1; the U log-weights are optimised by Nelder-Mead from several starts (bounded to
[-9, 9] by clipping).  No claim rests on this."""
import numpy as np
from scipy.optimize import minimize
from relax2 import Psi_matrix
from relax1 import frame


def lam_rel(logs, u):
    sU = np.exp(np.clip(logs, -9, 9))
    M = Psi_matrix(1.0, sU, 'V2')
    k = u + 1
    A = M[:-k, :-k]; B = M[:-k, -k:]; D = M[-k:, -k:]
    if np.linalg.eigvalsh(A).min() <= 0:
        return -1.0
    Gr = D - B.T @ np.linalg.solve(A, B)
    Pf = frame(1.0, sU)
    return np.linalg.eigvalsh(Pf.T @ Gr @ Pf).min()


if __name__ == '__main__':
    rng = np.random.default_rng(21)
    for u in (3, 4, 5):
        best = (9.0, None)
        for start in range(6):
            x0 = rng.uniform(-6, 6, u)
            r = minimize(lambda x: lam_rel(x, u), x0, method='Nelder-Mead',
                         options=dict(maxiter=250, xatol=1e-3, fatol=1e-6))
            if r.fun < best[0]:
                best = (r.fun, np.round(np.clip(r.x, -9, 9), 2))
        print(f'u = {u}: smallest lambda_min(G_rel) in the frame found: {best[0]:.5f} at log s_U = {best[1]} (sigma = 1)',
              flush=True)
