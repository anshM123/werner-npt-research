"""Section 11, F3: exact membership in K_1(3) by a dictionary.  Seed a dictionary with the terms of parametric fits
(6 + 12 and 12 + 40 filterings of rho + products) and 200 random columns; solve nonnegative least squares; then add
the filterings and products most aligned with the residual (local searches) and repeat.  An interior point of K_1(3)
is typically reproduced to machine precision, giving an explicit decomposition X = sum_k c_k col_k, c_k >= 0.
usage: python check_f3_nnls.py target.npy rounds   (target: a saved 3 x 3 operator; needs PyTorch and SciPy)
"""
import sys, numpy as np, torch
from scipy.optimize import nnls
torch.set_default_dtype(torch.float64)
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3): F[3 * i + j, 3 * j + i] = 1
rho = np.eye(9) - F / 2
Rt = torch.tensor(rho, dtype=torch.complex128)
X = np.load(sys.argv[1]); X = (X + X.conj().T) / 2; X = X / np.trace(X).real
ROUNDS = int(sys.argv[2])
def vecr(M): return np.concatenate([M.real.ravel(), M.imag.ravel()])
rng = np.random.default_rng(0)
cols = []
def add_filt(a, b):
    K = np.kron(a, b); Y = K @ rho @ K.conj().T; cols.append(Y / np.trace(Y).real)
def add_prod(x, y):
    v = np.kron(x, y); cols.append(np.outer(v, v.conj()) / np.vdot(v, v).real)
# seed with the terms of the best parametric fits (J filterings + P products)
def param_fit(J, P, RS):
    Yt = torch.tensor(X, dtype=torch.complex128); best = (np.inf, None)
    for r in range(RS):
        p = [(0.4 * torch.randn(*s, dtype=torch.complex128)).requires_grad_() for s in ((J, 3, 3), (J, 3, 3), (P, 3), (P, 3))]
        def model():
            V = torch.einsum('ka,kb->kab', p[2], p[3]).reshape(P, 9)
            Kt = torch.einsum('jab,jcd->jacbd', p[0], p[1]).reshape(J, 9, 9)
            return V.T @ V.conj() + (Kt @ Rt @ Kt.conj().transpose(1, 2)).sum(0)
        opt = torch.optim.LBFGS(p, lr=1, max_iter=8000, history_size=50, line_search_fn='strong_wolfe', tolerance_grad=1e-15, tolerance_change=1e-18)
        def closure():
            opt.zero_grad(); loss = torch.linalg.norm(model() - Yt) ** 2; loss.backward(); return loss
        opt.step(closure)
        with torch.no_grad():
            r_ = float(torch.linalg.norm(model() - Yt) / torch.linalg.norm(Yt))
            if r_ < best[0]: best = (r_, [q.detach().numpy().copy() for q in p])
    return best
torch.manual_seed(5)
for (J, P) in [(6, 12), (12, 40)]:
    r_, pp = param_fit(J, P, 6)
    print(f"seed fit J = {J}, P = {P}: {r_:.3e}", flush=True)
    for j in range(J): add_filt(pp[0][j], pp[1][j])
    for k in range(P): add_prod(pp[2][k], pp[3][k])
for _ in range(100):
    add_filt(rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3)), rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3)))
    add_prod(rng.standard_normal(3) + 1j * rng.standard_normal(3), rng.standard_normal(3) + 1j * rng.standard_normal(3))
def solve():
    Amat = np.stack([vecr(c) for c in cols], 1); coef, res = nnls(Amat, vecr(X), maxiter=20000)
    fitM = sum(c * M for c, M in zip(coef, cols))
    return coef, X - fitM
def best_columns(R, RS=16):
    """local searches maximising <R, col> over normalised filterings and products"""
    Rt_ = torch.tensor(R, dtype=torch.complex128); new = []
    for r in range(RS):
        a = torch.randn(3, 3, dtype=torch.complex128, requires_grad=True); b = torch.randn(3, 3, dtype=torch.complex128, requires_grad=True)
        def fv():
            K = torch.kron(a, b); Y = K @ Rt @ K.conj().T; return -torch.trace(Rt_ @ Y).real / torch.trace(Y).real
        opt = torch.optim.Adam([a, b], lr=0.03)
        for it in range(500):
            v = fv(); opt.zero_grad(); v.backward(); opt.step()
        new.append((-fv().item(), 'f', a.detach().numpy(), b.detach().numpy()))
        x = torch.randn(3, dtype=torch.complex128, requires_grad=True); y = torch.randn(3, dtype=torch.complex128, requires_grad=True)
        def gv():
            v = torch.kron(x, y); return -(v.conj() @ Rt_ @ v).real / (v.conj() @ v).real
        opt = torch.optim.Adam([x, y], lr=0.03)
        for it in range(300):
            v = gv(); opt.zero_grad(); v.backward(); opt.step()
        new.append((-gv().item(), 'p', x.detach().numpy(), y.detach().numpy()))
    return new
print(f"PT-min of X: {np.linalg.eigvalsh(X.reshape(3,3,3,3).transpose(0,3,2,1).reshape(9,9))[0]:+.4f}", flush=True)
for rd in range(ROUNDS):
    coef, R = solve()
    nR = np.linalg.norm(R) / np.linalg.norm(X)
    new = best_columns(R)
    gain = max(g for g, *_ in new)
    print(f"round {rd}: columns {len(cols)}, residual {nR:.3e}, best new alignment <R, col> = {gain:+.3e}", flush=True)
    if nR < 1e-12:
        print(f"exact decomposition: X is a nonnegative combination of {np.sum(coef > 0)} dictionary columns "
              f"(one-copy filterings and products) up to {nR:.1e}", flush=True); break
    if gain <= 1e-12:
        print("no improving column found: X appears to be outside K_1(3) at this distance", flush=True); break
    for g, kind, u, v in new:
        if g > 0: (add_filt if kind == 'f' else add_prod)(u, v)
