"""Minimum principle on a face: A_1 = w |psi><psi| (fixed unit psi, weight w free), A_2..A_n free, sum A = I.
Minimise the sector function mu_c over that face from random interior points of the face; report where runs end."""
import numpy as np, sys, warnings
warnings.filterwarnings("ignore")

from povm_form import sector_mu, m_povm
from scipy.optimize import minimize
n, c, starts, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]); rng = np.random.default_rng(seed)
psi = np.array([1.0, 0.4 + 0.3j]); psi /= np.linalg.norm(psi); P = np.outer(psi, psi.conj())
def face_povm(x):
    w = 1 / (1 + np.exp(-x[0]))                                   # weight of A_1 in (0, 1)
    G = (x[1:1 + 4 * (n - 1)] + 1j * x[1 + 4 * (n - 1):]).reshape(n - 1, 2, 2); B = [g @ g.conj().T for g in G]
    S = sum(B); rest = np.eye(2) - w * P                             # the other elements must sum to I - w P (positive definite for w < 1)
    ws, Vs = np.linalg.eigh(S); Si = Vs @ np.diag(ws ** -0.5) @ Vs.conj().T
    wr, Vr = np.linalg.eigh(rest); Rh = Vr @ np.diag(np.sqrt(wr)) @ Vr.conj().T
    return [w * P] + [Rh @ Si @ b @ Si @ Rh for b in B]
f = (lambda A: sector_mu(A, n, c)) if c >= 0 else (lambda A: m_povm(A, n))
for k in range(starts):
    x0 = rng.normal(size=1 + 8 * (n - 1))
    res = minimize(lambda x: f(face_povm(x)), x0, method='L-BFGS-B', options={'maxiter': 3000, 'gtol': 1e-12, 'ftol': 1e-15})
    A = face_povm(res.x); r1 = sum(1 for a in A[1:] if np.linalg.eigvalsh(a)[0] < 1e-6); w = np.trace(A[0]).real
    print(f"face n={n} c={c} start {k}: {f(face_povm(x0)):.5f} -> {res.fun:.6f}; weight of A_1 {w:.3f}; rank-<=1 among A_2..A_n {r1}/{n - 1}; "
          f"smallest eigenvalue among them {min(np.linalg.eigvalsh(a)[0] for a in A[1:]):.1e}", flush=True)
