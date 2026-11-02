"""Two-sector class in local dimension d (A_s perp B_s, dims <= 2): minimise q/||C||^2 (endpoint form, coefficient -1/2)
over rank two with ran C in A + B and arbitrary co-range.  usage: python check_other_dimensions.py d n trials seed"""
import sys, numpy as np, torch
torch.set_default_dtype(torch.float64)
d, n, TRIALS, SEED = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]); D = d ** n
torch.manual_seed(SEED); rng = np.random.default_rng(SEED); cplx = torch.complex128
wd = torch.eye(d, dtype=cplx).reshape(d * d) / np.sqrt(d)
def Phi(C):
    X = C.reshape((d,) * (2 * n)).permute([i for s in range(n) for i in (s, n + s)]).reshape((d * d,) * n)
    for s in range(n):
        tr = torch.tensordot(wd.conj(), X, dims=([0], [s])); X = X - (d / 2) * torch.moveaxis(torch.tensordot(wd, tr, dims=0), 0, s)
    return X.reshape([d, d] * n).permute([2 * s for s in range(n)] + [2 * s + 1 for s in range(n)]).reshape(D, D)
def q(C): return torch.vdot(C.reshape(-1), Phi(C).reshape(-1)).real
def kron_all(mats):
    out = torch.ones(1, 1, dtype=cplx)
    for M in mats: out = torch.kron(out, M)
    return out
worst = 1e9
for trial in range(TRIALS):
    dims = [[(2, 1), (1, 2), (1, 1), (2, 2)][rng.integers(4 if d >= 4 else 3)] for _ in range(n)]
    Us = []
    for s in range(n):
        Z = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)); Q_, _ = np.linalg.qr(Z); Us.append(torch.tensor(Q_))
    EA = kron_all([Us[s][:, :dims[s][0]] for s in range(n)]); EB = kron_all([Us[s][:, dims[s][0]:sum(dims[s])] for s in range(n)])
    E = torch.cat([EA, EB], 1); m = E.shape[1]
    def member(p):
        F = torch.complex(p[:2 * m], p[2 * m:4 * m]).reshape(m, 2); G = torch.complex(p[4 * m:4 * m + 2 * D], p[4 * m + 2 * D:]).reshape(D, 2)
        C = (E @ F) @ G.conj().T; return C / torch.linalg.norm(C)
    best = 1e9
    for st in range(3):
        p = torch.randn(4 * m + 4 * D, requires_grad=True)
        opt = torch.optim.LBFGS([p], lr=1, max_iter=2000, tolerance_grad=1e-14, tolerance_change=1e-18, history_size=50, line_search_fn='strong_wolfe')
        def closure():
            opt.zero_grad(); v = q(member(p)); v.backward(); return v
        opt.step(closure); best = min(best, q(member(p.detach())).item())
    worst = min(worst, best); print(f"trial {trial}: dims {dims}: min {best:+.2e}", flush=True)
print(f"d = {d}, n = {n}: min over trials {worst:+.2e}")
