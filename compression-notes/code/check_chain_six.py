"""Section 9: the next links of the 6 x 6 chain.  Are 6 x 6 local filterings of rho^{(x)M} (M = 3, 4) in
K_2(6) = cone{products, (a(x)b) rho^{(x)2} (a(x)b)^*, a, b : C^9 -> C^6}?  Targets: a two-copy NPT 6 x 6 filtering (gradient
path) padded with contractions on the extra copies, then moved by delta toward random M-copy filters.  Prints parametric
fits and saves the targets and fit terms (chain_target_M{M}_d{delta}.npy, chain_seeds_M{M}_d{delta}.pkl) for
check_chain_six_nnls.py.  usage: python check_chain_six.py M   (needs PyTorch and SciPy; about half an hour)
"""
import sys, itertools, numpy as np, torch
from scipy.optimize import nnls
torch.set_default_dtype(torch.float64)
M = int(sys.argv[1]); D = 6
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3): F[3 * i + j, 3 * j + i] = 1
rho = np.eye(9) - F / 2
R2 = np.kron(rho, rho).reshape([3] * 8).transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape(81, 81)
R2t = torch.tensor(R2, dtype=torch.complex128)
def filt(A, B, m):
    subsets = [S for r in range(m + 1) for S in itertools.combinations(range(m), r)]
    At, Bt = A.reshape([D] + [3] * m), B.reshape([D] + [3] * m); Y = 0
    for S in subsets:
        Sc = [s for s in range(m) if s not in S]
        pa = [0] + [1 + s for s in S] + [1 + s for s in Sc]
        Ap = At.permute(pa).reshape(D, 3 ** len(S), 3 ** len(Sc)); Bp = Bt.permute(pa).reshape(D, 3 ** len(S), 3 ** len(Sc))
        M1 = torch.einsum('axc,eyc->aexy', Ap, Ap.conj()); M2 = torch.einsum('byc,fxc->bfyx', Bp, Bp.conj())
        Y = Y + (-0.5) ** len(S) * torch.einsum('aexy,bfyx->abef', M1, M2).reshape(D * D, D * D)
    Y = (Y + Y.conj().T) / 2
    return Y / torch.trace(Y).real
def pt(Y): return Y.reshape(D, D, D, D).permute(0, 3, 2, 1).reshape(D * D, D * D)
# two-copy NPT seed filters (6 x 6), from a gradient path
torch.manual_seed(4)
A2 = torch.randn(D, 9, dtype=torch.complex128, requires_grad=True); B2 = torch.randn(D, 9, dtype=torch.complex128, requires_grad=True)
opt = torch.optim.Adam([A2, B2], lr=0.01)
for it in range(20000):
    Y = filt(A2, B2, 2); f = torch.linalg.eigvalsh(pt(Y))[0]
    if f.item() < -0.03: break
    opt.zero_grad(); f.backward(); opt.step()
print(f"two-copy seed: PT-min {f.item():+.4f}", flush=True)
u = torch.zeros(3 ** (M - 2), dtype=torch.complex128); u[0] = 1
A0 = torch.einsum('ax,y->axy', A2.detach(), u).reshape(D, 3 ** M); B0 = torch.einsum('ax,y->axy', B2.detach(), u).reshape(D, 3 ** M)
def param_fit(X, J, P, RS):
    Xt = torch.tensor(X, dtype=torch.complex128); best = (np.inf, None)
    for r in range(RS):
        p = [(0.3 * torch.randn(*s, dtype=torch.complex128)).requires_grad_() for s in ((J, D, 9), (J, D, 9), (P, D), (P, D))]
        def model():
            V = torch.einsum('ka,kb->kab', p[2], p[3]).reshape(P, D * D)
            Kt = torch.einsum('jab,jcd->jacbd', p[0], p[1]).reshape(J, D * D, 81)
            return V.T @ V.conj() + (Kt @ R2t @ Kt.conj().transpose(1, 2)).sum(0)
        o = torch.optim.LBFGS(p, lr=1, max_iter=6000, history_size=50, line_search_fn='strong_wolfe', tolerance_grad=1e-15, tolerance_change=1e-18)
        def cl():
            o.zero_grad(); loss = torch.linalg.norm(model() - Xt) ** 2; loss.backward(); return loss
        o.step(cl)
        with torch.no_grad():
            rr = float(torch.linalg.norm(model() - Xt) / torch.linalg.norm(Xt))
            if rr < best[0]: best = (rr, [q.detach().numpy().copy() for q in p])
    return best
def vecr(Mx): return np.concatenate([Mx.real.ravel(), Mx.imag.ravel()])
def nnls_check(X, seeds, extra=600):
    cols = []
    for pp in seeds:
        for j in range(pp[0].shape[0]):
            K = np.kron(pp[0][j], pp[1][j]); Y = K @ R2 @ K.conj().T; cols.append(Y / np.trace(Y).real)
        for k in range(pp[2].shape[0]):
            v = np.kron(pp[2][k], pp[3][k]); cols.append(np.outer(v, v.conj()) / np.vdot(v, v).real)
    g = np.random.default_rng(1)
    for _ in range(extra):
        a = g.standard_normal((D, 9)) + 1j * g.standard_normal((D, 9)); b = g.standard_normal((D, 9)) + 1j * g.standard_normal((D, 9))
        K = np.kron(a, b); Y = K @ R2 @ K.conj().T; cols.append(Y / np.trace(Y).real)
        x = g.standard_normal(D) + 1j * g.standard_normal(D); y = g.standard_normal(D) + 1j * g.standard_normal(D)
        v = np.kron(x, y); cols.append(np.outer(v, v.conj()) / np.vdot(v, v).real)
    Amat = np.stack([vecr(c) for c in cols], 1); coef, res = nnls(Amat, vecr(X), maxiter=50000)
    fitM = sum(c * Mc for c, Mc in zip(coef, cols))
    return np.linalg.norm(X - fitM) / np.linalg.norm(X), int(np.sum(coef > 0)), len(cols)
g2 = torch.Generator().manual_seed(9)
RA = torch.randn(D, 3 ** M, dtype=torch.complex128, generator=g2); RB = torch.randn(D, 3 ** M, dtype=torch.complex128, generator=g2)
torch.manual_seed(0)
for delta in [0.1, 0.25]:
    A = A0 / torch.linalg.norm(A0) + delta * RA / torch.linalg.norm(RA); B = B0 / torch.linalg.norm(B0) + delta * RB / torch.linalg.norm(RB)
    with torch.no_grad():
        X = filt(A, B, M).numpy()
    ptm = np.linalg.eigvalsh(pt(torch.tensor(X)).numpy())[0]
    f1 = param_fit(X, 16, 80, 3); f2 = param_fit(X, 32, 160, 2)
    rn, nact, ncols = nnls_check(X, [f1[1], f2[1]])
    np.save(f"chain_target_M{M}_d{delta}.npy", X)
    import pickle; pickle.dump([f1[1], f2[1]], open(f"chain_seeds_M{M}_d{delta}.pkl", "wb"))
    print(f"M = {M}, delta {delta}: PT-min {ptm:+.4f}; K_2(6) fits {f1[0]:.2e} (16+80), {f2[0]:.2e} (32+160); "
          f"dictionary NNLS {rn:.2e} ({nact} active of {ncols} columns)", flush=True)
