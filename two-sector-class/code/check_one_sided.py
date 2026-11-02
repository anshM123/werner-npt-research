"""One-sided two-sector class: check q(C) >= q(P_V C P_V) (V = A + B) for rank-two C with ran C in V and arbitrary
co-range, on random members and adversarially (minimise q(C) - q(P_V C P_V) relative to ||C||^2).
usage: python check_one_sided.py n trials seed"""
import sys, numpy as np, torch
torch.set_default_dtype(torch.float64)
n, TRIALS, SEED = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); D = 3 ** n; torch.manual_seed(SEED); rng = np.random.default_rng(SEED)
cplx = torch.complex128
w9 = torch.eye(3, dtype=cplx).reshape(9) / np.sqrt(3)
def Phi(C):
    X = C.reshape((3,) * (2 * n)).permute([i for s in range(n) for i in (s, n + s)]).reshape((9,) * n)
    for s in range(n):
        tr = torch.tensordot(w9.conj(), X, dims=([0], [s])); X = X - 1.5 * torch.moveaxis(torch.tensordot(w9, tr, dims=0), 0, s)
    return X.reshape([3, 3] * n).permute([2 * s for s in range(n)] + [2 * s + 1 for s in range(n)]).reshape(D, D)
def q(C): return torch.vdot(C.reshape(-1), Phi(C).reshape(-1)).real
def kron_all(mats):
    out = torch.ones(1, 1, dtype=cplx)
    for M in mats: out = torch.kron(out, M)
    return out
worst_rand, worst_adv = 1e9, 1e9
for trial in range(TRIALS):
    dims = [[(2, 1), (1, 2), (1, 1)][rng.integers(3)] for _ in range(n)]
    Us = []
    for s in range(n):
        Z = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); Q_, _ = np.linalg.qr(Z); Us.append(torch.tensor(Q_))
    EA = kron_all([Us[s][:, :dims[s][0]] for s in range(n)]); EB = kron_all([Us[s][:, dims[s][0]:sum(dims[s])] for s in range(n)])
    E = torch.cat([EA, EB], 1); m = E.shape[1]; PV = E @ E.conj().T
    def member(p):
        F = torch.complex(p[:2 * m], p[2 * m:4 * m]).reshape(m, 2); G = torch.complex(p[4 * m:4 * m + 2 * D], p[4 * m + 2 * D:]).reshape(D, 2)
        C = (E @ F) @ G.conj().T; return C / torch.linalg.norm(C)
    for k in range(10):
        C = member(torch.randn(4 * m + 4 * D)); worst_rand = min(worst_rand, (q(C) - q(PV @ C @ PV)).item())
    for st in range(3):
        p = torch.randn(4 * m + 4 * D, requires_grad=True)
        opt = torch.optim.LBFGS([p], lr=1, max_iter=2000, tolerance_grad=1e-14, tolerance_change=1e-18, history_size=50, line_search_fn='strong_wolfe')
        def closure():
            opt.zero_grad(); C = member(p); v = q(C) - q(PV @ C @ PV); v.backward(); return v
        opt.step(closure); C = member(p.detach()); worst_adv = min(worst_adv, (q(C) - q(PV @ C @ PV)).item())
print(f"n = {n}: min of q(C) - q(P_V C P_V) over random members {worst_rand:+.2e}, adversarial {worst_adv:+.2e}")
