"""Section 8: the smallest ratio q_n(C) / (theta(C) ||C||_F^2) over rank-two C = A B^* (floating point, Adam from random
complex factors).  theta(C) is the least eigenvalue of the one-copy marginals of the range and co-range projectors,
minimised over copies; it is zero exactly when C is thin at some copy.  usage: python thickness_margin.py n runs seed"""
import sys
import numpy as np
import torch

torch.set_default_dtype(torch.float64)
n, runs, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
d = 3 ** n


def q_seq(C):
    T = C.reshape([3] * (2 * n))
    X = T
    for s in range(n):
        tr = torch.diagonal(X, dim1=s, dim2=n + s).sum(-1)
        full = torch.einsum('ab,...->...ab', torch.eye(3, dtype=X.dtype), tr)
        X = X - 0.5 * torch.movedim(full, (2 * n - 2, 2 * n - 1), (s, n + s))
    return (T.conj() * X).sum().real


def proj(M):
    return M @ torch.linalg.solve(M.conj().T @ M, M.conj().T)


def lmin_marginal(P, s):
    others = [k for k in range(n) if k != s]
    Y = P.reshape([3] * (2 * n)).permute([s] + others + [n + s] + [n + k for k in others])
    rho = torch.einsum('aibi->ab', Y.reshape(3, 3 ** (n - 1), 3, 3 ** (n - 1)))
    return torch.linalg.eigvalsh((rho + rho.conj().T) / 2)[0]


def theta(A, B):
    PA, PB = proj(A), proj(B)
    return torch.stack([lmin_marginal(P, s) for P in (PA, PB) for s in range(n)]).min()


rng = np.random.default_rng(seed)
best = None
for r in range(runs):
    ps = [torch.tensor(rng.standard_normal((d, 2)), requires_grad=True) for _ in range(4)]
    opt = torch.optim.Adam(ps, lr=0.02)
    for it in range(6000):
        opt.zero_grad()
        A, B = torch.complex(ps[0], ps[1]), torch.complex(ps[2], ps[3])
        C = A @ B.conj().T
        ratio = q_seq(C) / (C.real ** 2 + C.imag ** 2).sum() / theta(A, B).clamp_min(1e-9)
        ratio.backward()
        opt.step()
    with torch.no_grad():
        A, B = torch.complex(ps[0], ps[1]), torch.complex(ps[2], ps[3])
        C = A @ B.conj().T
        qv = (q_seq(C) / (C.real ** 2 + C.imag ** 2).sum()).item()
        th = theta(A, B).item()
    print(f"run {r}: q/(theta |C|^2) = {qv / th:.4f}   q/|C|^2 = {qv:.3e}   theta = {th:.3e}", flush=True)
    best = qv / th if best is None else min(best, qv / th)
print(f"n = {n}: smallest ratio found {best:.4f}")
