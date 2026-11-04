"""Scope of the one-copy condition beyond the tangent class.
 (a) Bond-two class: range inside span{a_1..a_{s-1} (x) C^3 (x) b_{s+1}..b_n : s} (a = b is the tangent class;
     a = |0..0>, b = |1..1> contains the GHZ plane).  Balanced margins (scale-free) for random members, n = 3, 4, and
     near the GHZ plane.
 (b) Planes plus anchored single excitations, (x)_t B_t + sum_s C phi_s (x) f^(s) (contains the bond-two class): the stored
     operator data/anchored_pencil_psd_n3.npy is thick and has nu_s > 0 at every copy."""
import numpy as np, itertools
from scipy.optimize import minimize
from scipy.linalg import expm
from common import squeeze_forms, pencil_margin, Lmap, q, Bt, d
rng = np.random.default_rng(31)
B8 = Bt[[0, 1, 3, 4, 5, 6, 7, 8]]


def kronl(vs):
    v = vs[0]
    for w in vs[1:]: v = np.kron(v, w)
    return v


def act(A, X, n, s):
    D = d ** n; T = X.reshape([d] * (2 * n)); T = np.moveaxis(np.tensordot(A, T, axes=([1], [s])), 0, s)
    T = np.moveaxis(np.tensordot(T, A.conj().T, axes=([n + s], [0])), -1, n + s); return T.reshape(D, D)


def balanced_nu(X, n, s):
    def hm(h):
        H = np.einsum('k,kab->ab', h, B8); return H - np.trace(H) / 3 * np.eye(3)
    r = minimize(lambda h: np.log(np.linalg.norm(act(expm(hm(h)), X, n, s)) ** 2), np.zeros(8), method='BFGS', options={'gtol': 1e-9, 'maxiter': 2000})
    Y = act(expm(hm(r.x)), X, n, s); Y /= np.linalg.norm(Y); return pencil_margin(*squeeze_forms(Y, n, s))[0]


def bond_two_basis(a, b, n):
    cols = [kronl([a[t] for t in range(s)] + [np.eye(3)[k]] + [b[t] for t in range(s + 1, n)]) for s in range(n) for k in range(3)]
    Uq, S, _ = np.linalg.svd(np.array(cols).T, full_matrices=False); return Uq[:, :2 * n + 1]


cz = lambda *s: rng.normal(size=s) + 1j * rng.normal(size=s)
worst = -np.inf
print("(a) bond-two class, random members (balanced nu, minimum over copies)")
for n, samples in ((3, 6), (4, 3)):
    D = d ** n; vals = []
    for k in range(samples):
        a = [x / np.linalg.norm(x) for x in (cz(3) for _ in range(n))]; b = [x / np.linalg.norm(x) for x in (cz(3) for _ in range(n))]
        Bq = bond_two_basis(a, b, n); u = (Bq @ cz(2 * n + 1, 2)).T; v = cz(2, D)
        X = np.outer(u[0], v[0].conj()) + np.outer(u[1], v[1].conj()); X /= np.linalg.norm(X)
        vals.append(min(balanced_nu(X, n, s) for s in range(n)))
    worst = max(worst, max(vals))
    print(f"  n={n}: largest balanced min_s nu {max(vals):+.4f} over {samples} samples")
print("    near the GHZ plane span{|000>, |111>} (+ eps times a bond-two perturbation), n = 3")
n = 3; D = 27; e0, e1 = np.eye(3)[0], np.eye(3)[1]
for eps in (0.3, 0.1):
    pert = [sum(kronl([e0] * s + [cz(3)] + [e1] * (n - s - 1)) for s in range(n)) for _ in range(2)]
    u = np.array([kronl([e0] * n) + eps * pert[0], kronl([e1] * n) + eps * pert[1]]); v = cz(2, D)
    X = np.outer(u[0], v[0].conj()) + np.outer(u[1], v[1].conj()); X /= np.linalg.norm(X)
    nus = [balanced_nu(X, n, s) for s in range(n)]; worst = max(worst, max(nus))
    print(f"  eps {eps}: balanced nu per copy {np.round(nus, 4)}")
print("(b) planes plus anchored single excitations: stored operator")
X = np.load("../data/anchored_pencil_psd_n3.npy"); X /= np.linalg.norm(X)
Uu, S, _ = np.linalg.svd(X); u = (Uu[:, :2] * S[:2]).T
R = sum(np.outer(u[r], u[r].conj()) for r in range(2)).reshape([d] * 6); ratios = []
for s in range(3):
    T = np.moveaxis(R, [s, 3 + s], [0, 3]).reshape(3, 9, 3, 9); ev = np.linalg.eigvalsh(np.einsum('aibi->ab', T)); ratios.append(ev[0] / ev[-1])
nus = [pencil_margin(*squeeze_forms(X, 3, s))[0] for s in range(3)]
print(f"  local-support ratios {np.round(ratios, 3)}, nu per copy {np.round(nus, 4)}, R = {q(X, 3):.4f}")
assert worst < 0 and min(nus) > 0
print("RESULT: bond-two class nu < 0 in all tests; the anchored operator has nu > 0 at every copy (R > 0 there)")
