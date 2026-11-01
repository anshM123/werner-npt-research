"""Pair-traceless slice, better search: maximise S1/||C||^2 - M sum_{i<j} ||Tr_ij C||^2/||C||^2 with M ramped up,
from random starts and from perturbed product examples C = P_2 (x) |u><v| (u_j perp v_j), which give S1 = 2.
A value S1 > 2 with pair weight ~ 0 would be a counterexample to EP_n.  usage: python excl3.py n d starts"""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64)
n, d, NS = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); D = d ** n
def ptr(C, K):
    T = C.reshape([d] * (2 * n))
    for j in sorted(K, reverse=True):
        m = T.ndim // 2; T = torch.diagonal(T, dim1=j, dim2=j + m).sum(-1)
    return T
pairs = list(itertools.combinations(range(n), 2))
def prod_start(g):
    # A, B for C = P_2(site0) (x) |u><v|, u_j perp v_j, rotated randomly per site, plus noise
    A = torch.zeros(D, 2, dtype=torch.complex128); B = torch.zeros(D, 2, dtype=torch.complex128)
    for r in range(2):
        a = torch.zeros(d, dtype=torch.complex128); a[r] = 1
        va, vb = a, a.clone()
        for j in range(1, n):
            u = torch.zeros(d, dtype=torch.complex128); u[0] = 1; v = torch.zeros(d, dtype=torch.complex128); v[1] = 1
            va = torch.kron(va, u); vb = torch.kron(vb, v)
        A[:, r] = va; B[:, r] = vb
    A = A + 0.05 * torch.randn(D, 2, dtype=torch.complex128, generator=g); B = B + 0.05 * torch.randn(D, 2, dtype=torch.complex128, generator=g)
    return A, B
best = None
for s in range(NS):
    g = torch.Generator().manual_seed(500 + s)
    if s % 2 == 0: A, B = prod_start(g)
    else: A, B = torch.randn(D, 2, dtype=torch.complex128, generator=g), torch.randn(D, 2, dtype=torch.complex128, generator=g)
    A.requires_grad_(); B.requires_grad_()
    o = torch.optim.Adam([A, B], lr=0.01)
    for it in range(4000):
        M = 10.0 * (1 + it / 200.0)
        C = A @ torch.conj_physical(B).T; nn = torch.sum(torch.abs(C) ** 2)
        s1 = sum(torch.sum(torch.abs(ptr(C, [i])) ** 2) for i in range(n)) / nn
        s2 = sum(torch.sum(torch.abs(ptr(C, list(p))) ** 2) for p in pairs) / nn
        v = s1 - M * s2
        o.zero_grad(); (-v).backward(); o.step()
    print(f"start {s} ({'prod' if s % 2 == 0 else 'rand'}): S1 = {s1.item():.6f}, pair weight {s2.item():.2e}", flush=True)
    if best is None or s1.item() - 1e3 * s2.item() > best[0] - 1e3 * best[1]: best = (s1.item(), s2.item())
print(f"n={n} d={d}: best S1 = {best[0]:.6f} (pair weight {best[1]:.2e}); EP needs S1 <= 2 on the pair-traceless slice")
