"""Tomar Sec. 8 shape: range in span{0^k, 2_i 0^{k-1} (i in T0), 1^k, 2_j 1^{k-1} (j in T1)}, tails may share copies,
co-range ARBITRARY.  Adversarial min of q_n/|C|^2 over rank-2 C with that range.  usage: python check_tomar8_shared_tails.py n T0 T1 starts seed"""
import sys, numpy as np, torch
torch.set_default_dtype(torch.float64); cplx = torch.complex128
n = int(sys.argv[1]); T0 = [int(c) for c in sys.argv[2].split(',') if c]; T1 = [int(c) for c in sys.argv[3].split(',') if c]
ST, SEED = int(sys.argv[4]), int(sys.argv[5]); D = 3 ** n; torch.manual_seed(SEED)
w9 = torch.eye(3, dtype=cplx).reshape(9) / np.sqrt(3)
def Phi(C):
    X = C.reshape((3,) * (2 * n)).permute([i for s in range(n) for i in (s, n + s)]).reshape((9,) * n)
    for s in range(n):
        tr = torch.tensordot(w9.conj(), X, dims=([0], [s])); X = X - 1.5 * torch.moveaxis(torch.tensordot(w9, tr, dims=0), 0, s)
    return X.reshape([3, 3] * n).permute([2 * s for s in range(n)] + [2 * s + 1 for s in range(n)]).reshape(D, D)
q = lambda C: torch.vdot(C.reshape(-1), Phi(C).reshape(-1)).real
def basis_vec(digits):
    v = torch.zeros(D, dtype=cplx); v[np.ravel_multi_index(tuple(digits), (3,) * n)] = 1; return v
vecs = [basis_vec([0] * n), basis_vec([1] * n)]
for i in T0: d = [0] * n; d[i] = 2; vecs.append(basis_vec(d))
for j in T1: d = [1] * n; d[j] = 2; vecs.append(basis_vec(d))
Vb = torch.stack(vecs, 1); k = Vb.shape[1]
best = 1e9
for st in range(ST):
    p = torch.randn(4 * k + 4 * D, requires_grad=True)
    def member(p):
        A = Vb @ torch.complex(p[:2 * k], p[2 * k:4 * k]).reshape(k, 2)
        B = torch.complex(p[4 * k:4 * k + 2 * D], p[4 * k + 2 * D:]).reshape(D, 2)
        C = A @ B.conj().T; return C / torch.linalg.norm(C)
    opt = torch.optim.LBFGS([p], lr=1, max_iter=1500, tolerance_grad=1e-14, tolerance_change=1e-18, history_size=50, line_search_fn='strong_wolfe')
    def closure():
        opt.zero_grad(); v = q(member(p)); v.backward(); return v
    opt.step(closure)
    with torch.no_grad(): best = min(best, q(member(p)).item())
print(f"n={n} T0={T0} T1={T1}: adversarial min q_n/|C|^2 (any co-range) = {best:+.3e} over {ST} starts")
