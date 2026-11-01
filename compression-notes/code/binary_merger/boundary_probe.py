"""Boundary targets: (1) the 1-undistillability margin min over Schmidt-rank-2 unit psi of <psi|X^Gamma|psi>
(= 0 means X lies on the boundary of K_1(3), since K_1(3) members are 1-undistillable); (2) one-atom fits of
X_delta = X + delta I/9 (interior points of K_1(3) near X) for decreasing delta.  If SF failed at X, it would fail on a
neighbourhood (Sigma = {lam G + Sep} is closed), so X_delta would leave a robust residual for small delta.
usage: python boundary_probe.py file.npz index   (advmu_*.npz, or comm_s21_0.npz with index 0)"""
import sys, time, numpy as np, torch
from bm_core import G, RHO, ptmin, rand_c, D, N, pt
from bm_fit import fit_robust, norm_t
from ppt_margin import atom_hat_t
torch.set_num_threads(1)

fn, idx = sys.argv[1], int(sys.argv[2])
Z = np.load(fn, allow_pickle=True)
if 'X' in Z.files:
    X = norm_t(Z['X']); ins = [(Z['C'], Z['C']), (np.eye(D), np.eye(D))]
else:
    mode = str(Z['mode']); pp = Z['params'][idx]; tw = float(Z['tau'][idx]); w = np.exp(tw)
    if mode == 'nf':
        Gc = G(pp[0], pp[1]); X = norm_t(RHO / np.trace(RHO).real + w * Gc / np.trace(Gc).real); ins = [(pp[0], pp[1]), (np.eye(D), np.eye(D))]
    else:
        G1, G2 = G(pp[0], pp[1]), G(pp[2], pp[3]); X = norm_t(G1 / np.trace(G1).real + w * G2 / np.trace(G2).real); ins = [(pp[0], pp[1]), (pp[2], pp[3])]
g = np.random.default_rng(0)
XG = torch.tensor(pt(X), dtype=torch.complex128)
best = np.inf
for r in range(40):
    u = torch.randn(2, D, dtype=torch.complex128, requires_grad=True); v = torch.randn(2, D, dtype=torch.complex128, requires_grad=True)
    o = torch.optim.LBFGS([u, v], lr=1, max_iter=500, line_search_fn='strong_wolfe', tolerance_grad=1e-14, tolerance_change=1e-16)
    def cl():
        o.zero_grad(); psi = (torch.kron(u[0], v[0]) + torch.kron(u[1], v[1]))
        val = (psi.conj() @ XG @ psi).real / (psi.conj() @ psi).real; val.backward(); return val
    o.step(cl)
    with torch.no_grad():
        psi = (torch.kron(u[0], v[0]) + torch.kron(u[1], v[1])); best = min(best, float((psi.conj() @ XG @ psi).real / (psi.conj() @ psi).real))
print(f"{fn}[{idx}]: PT-min {ptmin(X):+.5f}, lambda_min {np.linalg.eigvalsh(X)[0]:.5f}; min over SR-2 psi of <psi|X^G|psi> = {best:+.3e}", flush=True)
for delta in (1e-2, 3e-3, 1e-3, 3e-4):
    t = time.time()
    Xd = norm_t(X + delta * np.eye(N) / N)
    out = fit_robust(Xd, P=80, restarts=3, rounds=30, seed=7, extra=ins)
    print(f"   delta {delta:.0e}: one-atom fit param {out['param']:.2e}, VP {out['vp']:.2e} (atom wt {out['weight']:.3f})  [{time.time() - t:.0f}s]", flush=True)
