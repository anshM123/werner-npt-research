"""Pair-traceless slice of the endpoint.  If Tr_ij C = 0 for all pairs i < j, then q_n(C) = ||C||^2 (1 - (n/2) a_1)
with a_1 = mean_i ||Tr_i C||^2/||C||^2, so EP_n on this slice <=> S_1 := sum_i ||Tr_i C||^2 <= 2 ||C||^2.
Maximise S_1/||C||^2 - M * (sum of pair weights) over rank-2 C; report S_1/||C||^2 and the exact q_n.
usage: python excl2.py n d M starts"""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64)
n, d, M, NS = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]); D = d ** n
def ptr(C, K):
    T = C.reshape([d] * (2 * n))
    for j in sorted(K, reverse=True):
        m = T.ndim // 2; T = torch.diagonal(T, dim1=j, dim2=j + m).sum(-1)
    return T
def qn(C):
    tot = 0
    for r in range(n + 1):
        for K in itertools.combinations(range(n), r):
            tot = tot + (-0.5) ** r * torch.sum(torch.abs(ptr(C, list(K))) ** 2)
    return tot
best = None
for s in range(NS):
    torch.manual_seed(100 + s)
    A = torch.randn(D, 2, dtype=torch.complex128, requires_grad=True); B = torch.randn(D, 2, dtype=torch.complex128, requires_grad=True)
    o = torch.optim.Adam([A, B], lr=0.02)
    for it in range(2500):
        C = A @ torch.conj_physical(B).T; nn = torch.sum(torch.abs(C) ** 2)
        s1 = sum(torch.sum(torch.abs(ptr(C, [i])) ** 2) for i in range(n)) / nn
        s2 = sum(torch.sum(torch.abs(ptr(C, [i, j])) ** 2) for i, j in itertools.combinations(range(n), 2)) / nn
        v = s1 - M * s2
        o.zero_grad(); (-v).backward(); o.step()
    with torch.no_grad():
        C = A @ torch.conj_physical(B).T; C = C / torch.linalg.norm(C)
        q = qn(C).item()
    print(f"start {s}: S1/||C||^2 = {s1.item():.5f}, pair weight {s2.item():.2e}, q_n/||C||^2 = {q:+.5f}", flush=True)
    if best is None or s1.item() > best[0]: best = (s1.item(), s2.item(), q)
print(f"n={n} d={d}: max S1 = {best[0]:.5f} (pair weight {best[1]:.2e}, q = {best[2]:+.5f}); EP on this slice needs S1 <= 2")
