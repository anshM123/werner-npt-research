"""Section 11, F3 at two copies (KL_1) from the dual side.  Fit T (the Tiles state, and Tiles with 5% noise) by K_1(3)
(10 filterings of rho + 30 products); with X* the fit, W = X* - T is, up to the fit error, a functional that is >= 0 on
K_1(3) and negative on T (Tiles is outside K_1(3): its range contains no product vector).  KL_1 says that such W stay
>= 0 on all 3 x 3 filterings of rho (x) rho.  Prints v_1(W) and v_2(W), the minima of <W, X>/Tr X over unit filters of one
and of two copies (local search, 24 restarts).  usage: python check_f3_dual.py   (needs PyTorch; about ten minutes)
"""
import numpy as np, torch
torch.set_default_dtype(torch.float64)
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3): F[3 * i + j, 3 * j + i] = 1
rho = np.eye(9) - F / 2
R1 = torch.tensor(rho, dtype=torch.complex128)
R2 = torch.tensor(np.kron(rho, rho).reshape([3] * 8).transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape(81, 81), dtype=torch.complex128)
e0, e1, e2 = np.eye(3)
upb = [np.kron(e0, (e0 - e1) / 2 ** .5), np.kron(e2, (e1 - e2) / 2 ** .5), np.kron((e0 - e1) / 2 ** .5, e2),
       np.kron((e1 - e2) / 2 ** .5, e0), np.kron((e0 + e1 + e2) / 3 ** .5, (e0 + e1 + e2) / 3 ** .5)]
Tiles = (np.eye(9) - sum(np.outer(x, x) for x in upb)) / 4
def project(T, J=10, P=30, RS=4):
    Tt = torch.tensor(T / np.trace(T).real, dtype=torch.complex128); best = (np.inf, None)
    for r in range(RS):
        p = [(0.4 * torch.randn(*s, dtype=torch.complex128)).requires_grad_() for s in ((J, 3, 3), (J, 3, 3), (P, 3), (P, 3))]
        def model():
            V = torch.einsum('ka,kb->kab', p[2], p[3]).reshape(P, 9)
            Kt = torch.einsum('jab,jcd->jacbd', p[0], p[1]).reshape(J, 9, 9)
            return V.T @ V.conj() + (Kt @ R1 @ Kt.conj().transpose(1, 2)).sum(0)
        opt = torch.optim.LBFGS(p, lr=1, max_iter=8000, history_size=50, line_search_fn='strong_wolfe', tolerance_grad=1e-15, tolerance_change=1e-18)
        def closure():
            opt.zero_grad(); loss = torch.linalg.norm(model() - Tt) ** 2; loss.backward(); return loss
        opt.step(closure)
        with torch.no_grad():
            M = model(); d = float(torch.linalg.norm(M - Tt))
            if d < best[0]: best = (d, (M - Tt).numpy())
    return best
def vmin(W, m, RS=24):
    Wt = torch.tensor(W, dtype=torch.complex128); R = R1 if m == 1 else R2; din = 3 ** m; best = np.inf
    for r in range(RS):
        A = torch.randn(3, din, dtype=torch.complex128, requires_grad=True); B = torch.randn(3, din, dtype=torch.complex128, requires_grad=True)
        def fv():
            K = torch.kron(A, B); X = K @ R @ K.conj().T; return torch.trace(Wt @ X).real / torch.trace(X).real
        opt = torch.optim.Adam([A, B], lr=0.03)
        for it in range(1500):
            v = fv(); opt.zero_grad(); v.backward(); opt.step()
        opt2 = torch.optim.LBFGS([A, B], lr=1, max_iter=300, line_search_fn='strong_wolfe', tolerance_grad=1e-14, tolerance_change=1e-16)
        def cl():
            opt2.zero_grad(); v = fv(); v.backward(); return v
        opt2.step(cl)
        best = min(best, fv().item())
    return best
torch.manual_seed(0)
for label, T in [("Tiles", Tiles), ("Tiles + 5% noise", 0.95 * Tiles + 0.05 * np.eye(9) / 9)]:
    d, W = project(T)
    W = (W + W.conj().T) / 2
    print(f"[{label}] distance to K_1(3) (fit): {d:.4f}; <W, T> = {np.trace(W @ T / np.trace(T)).real:+.4e}", flush=True)
    v1 = vmin(W, 1); v2 = vmin(W, 2)
    print(f"[{label}] v_1(W) = {v1:+.4e} (one-copy filterings), v_2(W) = {v2:+.4e} (two-copy filterings)", flush=True)
