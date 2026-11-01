"""Exclusive clicks: how large can the single-site partial-trace weight a_1 = mean_i ||Tr_i C||^2/||C||^2 be when
all pair weights ||Tr_ij C||^2 vanish, for rank-2 C on (C^d)^(x)n?  (Tomar's barrier profile is a = (1, 4/9, 0, ...).)
Maximise a_1 - M * a_2 for large M.  usage: python excl.py n d M starts"""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64)
n, d, M, NS = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]); D = d ** n
def ptr(C, K):
    T = C.reshape([d] * (2 * n))
    for j in sorted(K, reverse=True):
        m = T.ndim // 2; T = torch.diagonal(T, dim1=j, dim2=j + m).sum(-1)
    return T
best = None
for s in range(NS):
    torch.manual_seed(s)
    A = torch.randn(D, 2, dtype=torch.complex128, requires_grad=True); B = torch.randn(D, 2, dtype=torch.complex128, requires_grad=True)
    o = torch.optim.Adam([A, B], lr=0.02)
    for it in range(3000):
        C = A @ torch.conj_physical(B).T; nn = torch.sum(torch.abs(C) ** 2)
        a1 = sum(torch.sum(torch.abs(ptr(C, [i])) ** 2) for i in range(n)) / (n * nn)
        a2 = sum(torch.sum(torch.abs(ptr(C, [i, j])) ** 2) for i, j in itertools.combinations(range(n), 2)) / (n * (n - 1) / 2 * nn)
        v = a1 - M * a2
        o.zero_grad(); (-v).backward(); o.step()
    print(f"start {s}: a1 = {a1.item():.4f}, a2 = {a2.item():.2e}", flush=True)
    if best is None or v.item() > best[0]: best = (v.item(), a1.item(), a2.item())
print(f"n={n} d={d} M={M}: best a1 = {best[1]:.4f} with a2 = {best[2]:.2e}")
