"""Numerical evidence for Conjecture SC (single-copy squeeze condition) on the tangent class, and the contrast with the
two-excitation class.  nu_s = max_lam lambda_min(Q_s - lam N_s)/||N_s|| (lam in [-5, 5]); nu_s < 0 means the 9x9 pencil
at copy s is PSD for no lam, which is (stronger than) what the reduction needs.
 (a) random tangent operators, worst co-range for the range (n = 3, 4) or random co-range (n = 5);
 (b) targeted: project the known general single-copy minimum (orbit-notes/data/copywise_min_n3.npy) into the best
     tangent space and re-optimise inside the tangent class (maximise min_s nu_s): it runs to thin points, nu -> 0^-;
 (c) the same projection into the two-excitation class T^(2) gives a thick operator (data file) with nu_s > 0 at every copy."""
import numpy as np, itertools, sys
from scipy.optimize import minimize
from common import strings, Lmap, q, tangent_range, squeeze_forms, pencil_margin, d
rng = np.random.default_rng(11)


def worst_corange(u, n):
    D = d ** n; E = np.eye(D)
    cols = np.array([np.outer(u[r], E[j]) for r in range(2) for j in range(D)])
    imgs = np.array([Lmap(X, n) for X in cols]); Cm = cols.reshape(len(cols), -1); Im = imgs.reshape(len(cols), -1)
    Qm = Cm.conj() @ Im.T; Nm = Cm.conj() @ Cm.T; Qm = (Qm + Qm.conj().T) / 2; Nm = (Nm + Nm.conj().T) / 2
    w, P = np.linalg.eigh(Nm); keep = w > 1e-10 * w.max(); Pk = P[:, keep] / np.sqrt(w[keep])
    ev, evec = np.linalg.eigh(Pk.conj().T @ Qm @ Pk); c = Pk @ evec[:, 0]
    X = np.tensordot(c, cols, axes=(0, 0)); return X / np.linalg.norm(X), ev[0]


def local_ratios(u, n):
    R = sum(np.outer(u[r], u[r].conj()) for r in range(u.shape[0])).reshape([d] * (2 * n)); out = []
    for s in range(n):
        T = np.moveaxis(R, [s, n + s], [0, n]).reshape(d, d ** (n - 1), d, d ** (n - 1))
        ev = np.linalg.eigvalsh(np.einsum('aibi->ab', T)); out.append(ev[0] / ev[-1])
    return np.array(out)


print("(a) random tangent operators")
worst_a = -np.inf
for n, trials, worst in ((3, 6, True), (4, 3, True), (5, 3, False)):
    D = d ** n
    for t in range(trials):
        U = [rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)) for _ in range(n)]; a = rng.normal(size=2) + 1j * rng.normal(size=2)
        u = tangent_range(U, a, n)
        if worst:
            X, R = worst_corange(u, n)
        else:
            v = rng.normal(size=(2, D)) + 1j * rng.normal(size=(2, D)); X = np.outer(u[0], v[0].conj()) + np.outer(u[1], v[1].conj())
            X /= np.linalg.norm(X); R = q(X, n)
        nus = [pencil_margin(*squeeze_forms(X, n, s))[0] for s in range(n)]
        worst_a = max(worst_a, max(nus))
        print(f"  n={n} trial {t}: R = {R:.4f}, nu per copy {np.round(nus, 4)}, local-support ratios {np.round(local_ratios(u, n), 3)}")
print(f"  largest nu over (a): {worst_a:+.4f}")

print("(b) targeted test inside the tangent class, n = 3")
n = 3; D = d ** n; st = strings(n)
C0 = np.load("../../orbit-notes/data/copywise_min_n3.npy").reshape(D, D)
Uc, S, Vh = np.linalg.svd(C0); P = Uc[:, :2] @ Uc[:, :2].conj().T
allowed = np.array([sum(c != 0 for c in s_) <= 1 for s_ in st]); na = allowed.sum()
def lu(x):
    Us = [np.linalg.qr(x[18 * s:18 * s + 9].reshape(3, 3) + 1j * x[18 * s + 9:18 * s + 18].reshape(3, 3))[0] for s in range(n)]
    W = Us[0]
    for t in range(1, n): W = np.kron(W, Us[t])
    return W
best = (0, None)
for k in range(10):
    r = minimize(lambda x: -np.trace((lu(x).conj().T @ P @ lu(x))[np.ix_(allowed, allowed)]).real / 2, rng.normal(size=18 * n), method='BFGS', options={'maxiter': 300})
    if -r.fun > best[0]: best = (-r.fun, r.x)
W = lu(best[1]); A = W.conj().T @ Uc[:, :2] * S[:2]; Bv = W.conj().T @ Vh[:2].conj().T; A[~allowed] = 0
print(f"  fraction of the range of the general single-copy minimum inside the best tangent space: {best[0]:.4f}")
def build(x):
    a_ = (x[:2 * na] + 1j * x[2 * na:4 * na]).reshape(na, 2); b_ = (x[4 * na:4 * na + 2 * D] + 1j * x[4 * na + 2 * D:]).reshape(D, 2)
    Au = np.zeros((D, 2), complex); Au[allowed] = a_; X = Au @ b_.conj().T; return X / np.linalg.norm(X), Au.T
def vals(x):
    X, u = build(x); return np.array([pencil_margin(*squeeze_forms(X, n, s))[0] for s in range(n)]), u
x0 = np.concatenate([A[allowed].real.ravel(), A[allowed].imag.ravel(), Bv.real.ravel(), Bv.imag.ravel()])
def obj(x):
    m, _ = vals(x); return -(m.min() - np.log(np.sum(np.exp(-200 * (m - m.min())))) / 200)
v0, u0 = vals(x0); print(f"  after projection: nu {np.round(v0, 4)}, local-support ratios {np.round(local_ratios(u0, n), 4)}")
r = minimize(obj, x0, method='L-BFGS-B', options={'maxiter': 300, 'maxfun': 40000})
v1, u1 = vals(r.x); print(f"  after re-optimisation inside the tangent class: nu {np.round(v1, 6)}, local-support ratios {np.round(local_ratios(u1, n), 5)}")

print("(c) two-excitation class: the stored operator")
X2 = np.load("../data/two_excitation_pencil_psd_n3.npy"); n = 3
Uq, S2, _ = np.linalg.svd(X2); u2 = (Uq[:, :2] * S2[:2]).T
outside = np.array([sum(c != 0 for c in s_) > 2 for s_ in st])
print(f"  rank {np.sum(S2 > 1e-9 * S2[0])}, range weight outside T^(2): {np.sum(abs(Uq[outside, :2]) ** 2):.1e}, "
      f"local-support ratios {np.round(local_ratios(u2, n), 3)}")
nus2 = [pencil_margin(*squeeze_forms(X2, n, s))[0] for s in range(n)]
print(f"  nu per copy {np.round(nus2, 4)} (> 0: the pencil is positive definite for some lam at every copy)")
assert worst_a < 0 and max(nus2) > 0 and min(nus2) > 0
print("RESULT: tangent class nu < 0 in all tests; two-excitation operator has nu > 0 at every copy")
