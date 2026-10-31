"""Section 9: exact membership in K_2(D) by column generation.  Dictionary: the fit terms of check_chain_six.py, 600 random
filterings of rho^{(x)2} and 600 products; nonnegative least squares; then add the filterings and products most aligned
with the residual (24 local searches each per round).  usage: python check_chain_six_nnls.py chain_target_M4.npy
chain_seeds_M4.pkl 60   (needs PyTorch and SciPy; about forty minutes)
"""
import sys, pickle, numpy as np, torch
from scipy.optimize import nnls
torch.set_default_dtype(torch.float64)
X = np.load(sys.argv[1]); X = (X + X.conj().T) / 2; X = X / np.trace(X).real; D = int(round(X.shape[0] ** 0.5))
seeds = pickle.load(open(sys.argv[2], 'rb')); ROUNDS = int(sys.argv[3])
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3): F[3 * i + j, 3 * j + i] = 1
rho = np.eye(9) - F / 2
R2 = np.kron(rho, rho).reshape([3] * 8).transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape(81, 81); R2t = torch.tensor(R2, dtype=torch.complex128)
cols = []
def add_f(a, b):
    K = np.kron(a, b); Y = K @ R2 @ K.conj().T; cols.append(Y / np.trace(Y).real)
def add_p(x, y):
    v = np.kron(x, y); cols.append(np.outer(v, v.conj()) / np.vdot(v, v).real)
for pp in seeds:
    for j in range(pp[0].shape[0]): add_f(pp[0][j], pp[1][j])
    for k in range(pp[2].shape[0]): add_p(pp[2][k], pp[3][k])
g = np.random.default_rng(1)
for _ in range(600):
    add_f(g.standard_normal((D, 9)) + 1j * g.standard_normal((D, 9)), g.standard_normal((D, 9)) + 1j * g.standard_normal((D, 9)))
    add_p(g.standard_normal(D) + 1j * g.standard_normal(D), g.standard_normal(D) + 1j * g.standard_normal(D))
def vecr(Mx): return np.concatenate([Mx.real.ravel(), Mx.imag.ravel()])
for rd in range(ROUNDS):
    Amat = np.stack([vecr(c) for c in cols], 1); coef, _ = nnls(Amat, vecr(X), maxiter=100000)
    R = X - sum(c * Mc for c, Mc in zip(coef, cols)); nR = np.linalg.norm(R) / np.linalg.norm(X)
    Rt_ = torch.tensor(R, dtype=torch.complex128); new = []
    for r in range(24):
        a = torch.randn(D, 9, dtype=torch.complex128, requires_grad=True); b = torch.randn(D, 9, dtype=torch.complex128, requires_grad=True)
        o = torch.optim.Adam([a, b], lr=0.03)
        for it in range(400):
            K = torch.kron(a, b); Y = K @ R2t @ K.conj().T; v = -torch.trace(Rt_ @ Y).real / torch.trace(Y).real
            o.zero_grad(); v.backward(); o.step()
        new.append((-v.item(), 'f', a.detach().numpy(), b.detach().numpy()))
        x = torch.randn(D, dtype=torch.complex128, requires_grad=True); y = torch.randn(D, dtype=torch.complex128, requires_grad=True)
        o = torch.optim.Adam([x, y], lr=0.03)
        for it in range(300):
            w = torch.kron(x, y); v = -(w.conj() @ Rt_ @ w).real / (w.conj() @ w).real
            o.zero_grad(); v.backward(); o.step()
        new.append((-v.item(), 'p', x.detach().numpy(), y.detach().numpy()))
    gain = max(gg for gg, *_ in new)
    print(f"round {rd}: columns {len(cols)}, active {int(np.sum(coef > 0))}, residual {nR:.3e}, best alignment {gain:+.3e}", flush=True)
    if nR < 1e-12: print("exact decomposition", flush=True); break
    if gain <= 1e-13: print(f"no improving column: residual {nR:.3e} estimates the distance from K_2({D})", flush=True); break
    for gg, kind, u, v in new:
        if gg > 0: (add_f if kind == 'f' else add_p)(u, v)
