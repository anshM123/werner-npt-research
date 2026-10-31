"""Structure of KL_1 decompositions: for NPT 3x3 filterings X of rho (x) rho, fit X = (a(x)b) rho (a(x)b)^* + products
   (ONE filtering), then ask whether the filtering is a contraction of the target's filters: is a a^* proportional to
   A_u A_u^* for some u, with A_u = A(I (x) |u>) (copy-2 contraction) or A(|u> (x) I) (copy-1 contraction)?
   (a is only determined up to a -> a U with U unitary, since rho is U (x) U invariant, so compare a a^*.)"""
import numpy as np, torch
torch.set_default_dtype(torch.float64)
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3): F[3 * i + j, 3 * j + i] = 1
rho = np.eye(9) - F / 2
rho2 = np.kron(rho, rho).reshape([3] * 8).transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape(81, 81)
R1 = torch.tensor(rho, dtype=torch.complex128); R2 = torch.tensor(rho2, dtype=torch.complex128)
def pt(Y): return Y.reshape(3, 3, 3, 3).permute(0, 3, 2, 1).reshape(9, 9)
def fit1(Y, P=12, RS=10, IT=4000):
    Yt = torch.tensor(Y, dtype=torch.complex128); best = (np.inf, None)
    for _ in range(RS):
        p = [(0.5 * torch.randn(*s, dtype=torch.complex128)).requires_grad_() for s in ((3, 3), (3, 3), (P, 3), (P, 3))]
        def model():
            V = torch.einsum('ka,kb->kab', p[2], p[3]).reshape(P, 9)
            K = torch.kron(p[0], p[1]); return V.T @ V.conj() + K @ R1 @ K.conj().T
        opt = torch.optim.LBFGS(p, lr=1, max_iter=IT, history_size=50, line_search_fn='strong_wolfe', tolerance_grad=1e-15, tolerance_change=1e-18)
        def closure():
            opt.zero_grad(); loss = torch.linalg.norm(model() - Yt) ** 2; loss.backward(); return loss
        opt.step(closure)
        with torch.no_grad():
            r = float(torch.linalg.norm(model() - Yt) / torch.linalg.norm(Yt))
            if r < best[0]: best = (r, (p[0].detach().numpy().copy(), p[1].detach().numpy().copy()))
    return best
def contraction_gap(a, A, which):
    """min over unit u of || a a^*/Tr - A_u A_u^*/Tr ||, A_u = contraction of A (3x9) on copy 'which' with u."""
    M = a @ a.conj().T; M = M / np.trace(M).real
    At = A.reshape(3, 3, 3)
    best = np.inf
    for s in range(40):
        u = np.random.standard_normal(3) + 1j * np.random.standard_normal(3)
        ut = torch.tensor(u, requires_grad=True)
        Att = torch.tensor(At)
        opt = torch.optim.Adam([ut], lr=0.05)
        for it in range(400):
            uu = ut / torch.linalg.norm(ut)
            Au = torch.einsum('oij,j->oi', Att, uu) if which == 2 else torch.einsum('oij,i->oj', Att, uu)
            N = Au @ Au.conj().T; N = N / torch.trace(N).real
            g = torch.linalg.norm(N - torch.tensor(M))
            opt.zero_grad(); g.backward(); opt.step()
        best = min(best, g.item())
    return best
torch.manual_seed(7); np.random.seed(7)
for seed in range(4):
    A = torch.randn(3, 9, dtype=torch.complex128, requires_grad=True); B = torch.randn(3, 9, dtype=torch.complex128, requires_grad=True)
    opt = torch.optim.Adam([A, B], lr=0.01)
    for it in range(20000):
        K = torch.kron(A, B); Y = K @ R2 @ K.conj().T; Y = (Y + Y.conj().T) / 2; Y = Y / torch.trace(Y).real
        f = torch.linalg.eigvalsh(pt(Y))[0]
        if f.item() < -0.01: break
        opt.zero_grad(); f.backward(); opt.step()
    Yn = Y.detach().numpy(); An, Bn = A.detach().numpy(), B.detach().numpy()
    r, (a, b) = fit1(Yn)
    ga2, ga1 = contraction_gap(a, An, 2), contraction_gap(a, An, 1)
    gb2, gb1 = contraction_gap(b, Bn, 2), contraction_gap(b, Bn, 1)
    sa = np.linalg.svd(a, compute_uv=False); sA = np.linalg.svd(An, compute_uv=False)
    print(f"target {seed} (PT-min {f.item():+.4f}): one-filtering fit {r:.2e} | a vs copy-2/copy-1 contractions of A: "
          f"{ga2:.2e} / {ga1:.2e} | b: {gb2:.2e} / {gb1:.2e} | sv(a)/max {np.round(sa / sa[0], 3)}", flush=True)
