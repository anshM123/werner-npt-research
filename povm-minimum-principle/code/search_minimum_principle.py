"""Numerical tests of the minimum principle.
  descent S n starts seed      : minimise m over n-outcome qubit POVMs (single-excitation range) from random interior points
  descent T n starts seed      : minimise m over tangent planes (with a vacuum amplitude) from random points
  descent c n starts seed C    : minimise the sector function mu_C over POVMs
  adversary S n starts seed h delta [C] : maximise phi_h(A) = min over directions Delta of [f(A + h Delta) - f(A)]/h over
                                 POVMs with every A_s >= delta (f = m, or mu_C); phi_h > 0 would be a strict interior
                                 local minimum at scale h.  Directions: +/- an orthonormal basis of the tangent space
                                 {sum Delta_s = 0} and 40 random unit directions.
Each descent reports where it ends: the smallest eigenvalue among the A_s (-> 0 at the boundary) and how many A_s
have rank <= 1 (all of them: an extreme point)."""
import numpy as np, sys, warnings
from scipy.optimize import minimize
from povm_form import m_povm, sector_mu, povm_from, tangent_from, m_tangent
warnings.filterwarnings("ignore")
kind, sub, n, starts, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
rng = np.random.default_rng(seed); minEig = lambda A: min(np.linalg.eigvalsh(a)[0] for a in A)
if kind == "descent":
    C = int(sys.argv[6]) if sub == "c" else None
    for k in range(starts):
        if sub == "T":
            x0 = rng.normal(size=2 * (4 * n + 2)); f = lambda x: m_tangent(x, n)
        else:
            x0 = rng.normal(size=8 * n); f = (lambda x: sector_mu(povm_from(x, n), n, C)) if sub == "c" else (lambda x: m_povm(povm_from(x, n), n))
        res = minimize(f, x0, method='L-BFGS-B', options={'maxiter': 4000, 'gtol': 1e-12, 'ftol': 1e-15})
        if sub == "T":
            U, a = tangent_from(res.x, n); A = [u.conj().T @ u for u in U]; extra = f", vacuum weight {np.vdot(a, a).real:.3f}"
        else:
            A = povm_from(res.x, n); extra = ""
        r1 = sum(1 for a_ in A if np.linalg.eigvalsh(a_)[0] < 1e-6)
        print(f"descent {sub} n={n}{'' if C is None else f' c={C}'} start {k}: {f(x0):.5f} -> {res.fun:.6f}; smallest eigenvalue among A_s {minEig(A):.1e}; "
              f"rank-<=1 elements {r1}/{n}{extra}", flush=True)
else:
    h, delta = float(sys.argv[6]), float(sys.argv[7]); C = int(sys.argv[8]) if len(sys.argv) > 8 else None
    f = (lambda A: sector_mu(A, n, C)) if C is not None else (lambda A: m_povm(A, n))
    herm = [np.array([[1, 0], [0, 0]]), np.array([[0, 0], [0, 1]]), np.array([[0, 1], [1, 0]]) / np.sqrt(2), np.array([[0, -1j], [1j, 0]]) / np.sqrt(2)]
    raw = []
    for s in range(n - 1):
        for hk in herm:
            Dl = [np.zeros((2, 2), complex) for _ in range(n)]; Dl[s] = hk.astype(complex); Dl[n - 1] = -hk.astype(complex)
            raw.append(np.concatenate([np.concatenate([d.real.ravel(), d.imag.ravel()]) for d in Dl]))
    Q, _ = np.linalg.qr(np.array(raw).T); T = Q.shape[1]
    vecs = [Q[:, k] for k in range(T)] + [-Q[:, k] for k in range(T)] + [Q @ (lambda v: v / np.linalg.norm(v))(rng.normal(size=T)) for _ in range(40)]
    dirs = [[(lambda d: (d + d.conj().T) / 2)((v[8 * s:8 * s + 4] + 1j * v[8 * s + 4:8 * s + 8]).reshape(2, 2)) for s in range(n)] for v in vecs]
    phi = lambda A: (lambda m0: min((f([a + h * d for a, d in zip(A, Dl)]) - m0) / h for Dl in dirs))(f(A))
    best = -np.inf
    for k in range(starts):
        obj = lambda x: (lambda A: -phi(A) + 1e3 * max(0.0, delta - minEig(A)) ** 2)(povm_from(x, n))
        res = minimize(obj, rng.normal(size=8 * n), method='Nelder-Mead', options={'maxiter': 6000, 'xatol': 1e-6, 'fatol': 1e-10})
        A = povm_from(res.x, n); v = phi(A)
        if minEig(A) > 0.9 * delta: best = max(best, v)
        print(f"adversary n={n}{'' if C is None else f' c={C}'} h={h} start {k}: max phi_h {v:+.3e} at value {f(A):.5f}, min eig A_s {minEig(A):.3f}", flush=True)
    print(f"adversary n={n}: largest phi_h at interior points {best:+.3e} (> 0 would be a strict interior local minimum)")
