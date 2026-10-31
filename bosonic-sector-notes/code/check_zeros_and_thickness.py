"""Section 3, continued (local optimisation; not certified).
   (a) zeros at n = 3: singular values, and the spectrum of the one-site marginal of P_ran + P_coran;
   (b) thickness-constrained minima at n = 2, 3, and the best qubit compression C_z = E_0 C E_0 of each minimiser."""
import numpy as np
from scipy.optimize import minimize
from bosonic import basis, ann, q_and_grad

rng = np.random.default_rng(5)

def unpack(x, D):
    Y = (x[:2 * D] + 1j * x[2 * D:4 * D]).reshape(D, 2); Z = (x[4 * D:6 * D] + 1j * x[6 * D:]).reshape(D, 2)
    return Y, Z

def marg(M, n):
    A = ann(n)
    return np.array([[np.trace(A[i] @ M @ A[j].T) / n for j in range(3)] for i in range(3)])

def R(C, n, t=0.5):
    return q_and_grad(C, n, t)[0] / np.vdot(C, C).real

def planes(C):
    U, s, Vh = np.linalg.svd(C); u = U[:, :2]; v = Vh[:2].conj().T
    return s, u @ u.conj().T, v @ v.conj().T

print("(a) zeros at n = 3")
n = 3; D = len(basis(n))
def f(x):
    Y, Z = unpack(x, D); C = Y @ Z.conj().T; nn = np.vdot(C, C).real; q, G = q_and_grad(C, n)
    val = q / nn; Gc = (G - val * C) / nn; GY = Gc @ Z; GZ = Gc.conj().T @ Y
    return val, 2 * np.concatenate([GY.real.ravel(), GY.imag.ravel(), GZ.real.ravel(), GZ.imag.ravel()])
for k in range(8):
    r = minimize(f, rng.standard_normal(8 * D), jac=True, method='L-BFGS-B',
                 options={'maxiter': 50000, 'gtol': 1e-13, 'ftol': 1e-16})
    Y, Z = unpack(r.x, D); C = Y @ Z.conj().T; C /= np.linalg.norm(C)
    s, Pu, Pv = planes(C); ev = np.linalg.eigvalsh(marg(Pu + Pv, n))[::-1]
    print(f"  q/||C||^2 = {r.fun:+.1e}   singular values {np.round(s[:2], 4)}   "
          f"marginal of P_ran + P_coran (normalised) {np.round(ev / ev.sum(), 4)}")

print("(b) thickness-constrained minima and their best qubit compressions")
def Ez(z, n):
    A = ann(n); z = z / np.linalg.norm(z)
    N = sum(z[i] * np.conj(z[j]) * (A[i].T @ A[j]) for i in range(3) for j in range(3))
    w, V = np.linalg.eigh((N + N.conj().T) / 2); P = V[:, w < 0.5]; return P @ P.conj().T
def thick(C, n):
    s, Pu, Pv = planes(C); return np.linalg.eigvalsh(marg(Pu + Pv, n) / 4)[0]
def best_z(C, n):
    def val(x):
        E = Ez(x[:3] + 1j * x[3:], n); Cz = E @ C @ E; nz = np.vdot(Cz, Cz).real
        return q_and_grad(Cz, n)[0] / nz if nz > 1e-10 else 1e3
    X = rng.standard_normal((300, 6)); vals = np.array([val(x) for x in X])
    return min(minimize(val, X[i], method='Nelder-Mead', options={'maxiter': 1500, 'xatol': 1e-9, 'fatol': 1e-12}).fun
               for i in np.argsort(vals)[:5])
for n, taus in ((2, (0.01, 0.06, 0.15)), (3, (0.01, 0.06, 0.15))):
    D = len(basis(n))
    for tau in taus:
        g = lambda x: (lambda C: R(C, n) + 200.0 * max(0.0, tau - thick(C, n)) ** 2)(
            (lambda Y, Z: Y @ Z.conj().T)(*unpack(x, D)))
        r = min((minimize(g, rng.standard_normal(8 * D), method='L-BFGS-B', options={'maxiter': 3000})
                 for _ in range(2)), key=lambda r: r.fun)
        Y, Z = unpack(r.x, D); C = Y @ Z.conj().T
        rc, th, bz = R(C, n), thick(C, n), best_z(C, n)
        print(f"  n={n} tau={tau:.2f}: thickness {th:.4f}  R(C) {rc:+.5f}  min_z R(C_z) {bz:+.5f}  gap {bz - rc:+.2e}",
              flush=True)
