"""ORB_m test for GENERAL range planes: m(U) = min over co-ranges of R = lambda_min of K in the Gram-whitened frame,
minimised over the range-side filter orbit U -> (x)_s exp(H_s) U (H_s Hermitian 3x3; unitaries act trivially).
A run ends at the boundary if some copy becomes thin (range local support ratio -> 0).  usage: python orbm_general.py n starts seed"""
import numpy as np, sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "../../povm-minimum-principle/code")
from povm_form import psi_n
from scipy.linalg import expm
from scipy.optimize import minimize
n, starts, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); rng = np.random.default_rng(seed); D = 3 ** n
gm = []
for i in range(3):
    E = np.zeros((3, 3), complex); E[i, i] = 1; gm.append(E)
for i in range(3):
    for j in range(i + 1, 3):
        E = np.zeros((3, 3), complex); E[i, j] = E[j, i] = 1; gm.append(E)
        E = np.zeros((3, 3), complex); E[i, j] = -1j; E[j, i] = 1j; gm.append(E)
gm = np.array(gm)
def m_plane(u):
    if not np.all(np.isfinite(u)): return 10.0
    u = u / np.linalg.norm(u)
    G = np.array([[np.vdot(u[r], u[rp]) for rp in range(2)] for r in range(2)])
    if np.linalg.cond(G) > 1e10: return 10.0
    L = np.linalg.cholesky(G); Li = np.linalg.inv(L)
    w = Li @ u                                                   # orthonormal basis of the plane
    K = np.zeros((2 * D, 2 * D), complex)
    for r in range(2):
        for rp in range(2):
            K[r * D:(r + 1) * D, rp * D:(rp + 1) * D] = psi_n(np.outer(w[r], w[rp].conj()), n)
    try:
        return np.linalg.eigvalsh((K + K.conj().T) / 2)[0]
    except np.linalg.LinAlgError:
        return 10.0
def filt(x):
    G = np.eye(1)
    for s in range(n):
        G = np.kron(G, expm(np.einsum('k,kab->ab', x[9 * s:9 * s + 9], gm)))
    return G
def thin_ratio(u):
    P = sum(np.outer(v, v.conj()) for v in u); out = []
    for s in range(n):
        T = np.moveaxis(P.reshape([3] * (2 * n)), [s, n + s], [0, n]).reshape(3, D // 3, 3, D // 3)
        e = np.linalg.eigvalsh(np.einsum('aibi->ab', T)); out.append(e[0] / e[-1])
    return min(out)
for k in range(starts):
    u0 = rng.normal(size=(2, D)) + 1j * rng.normal(size=(2, D))
    f = lambda x: m_plane(np.array([filt(x) @ v for v in u0]))
    res = minimize(f, np.zeros(9 * n), method='L-BFGS-B', bounds=[(-8, 8)] * (9 * n), options={'maxiter': 3000, 'gtol': 1e-11, 'ftol': 1e-14})
    uf = np.array([filt(res.x) @ v for v in u0])
    print(f"ORB_m general n={n} start {k}: m {f(np.zeros(9 * n)):.5f} -> {res.fun:.6f}; smallest range local-support ratio {thin_ratio(uf):.1e}; "
          f"max |x| {np.abs(res.x).max():.1f}", flush=True)
