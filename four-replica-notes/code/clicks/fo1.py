"""FO1: for rank-<=2 R on (C^d)^(x)n, sum_i W_{[n]-i}(R) <= (2/d) ||R||_F^2, where W_{[n]-i} is the weight of R in the
sector 'identity at site i, traceless at every other site'.  W_{[n]-i}(R) = ||Pi_tl(Tr_i R)||^2 / d, and for an
operator X on m sites ||Pi_tl X||^2 = sum_K (-1/d)^{|K|} ||Tr_K X||^2.  Adversarial maximisation of the ratio.
usage: python fo1.py n d starts"""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64)
n, d, NS = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); D = d ** n
def ptr(C, K):
    T = C.reshape([d] * (2 * n))
    for j in sorted(K, reverse=True):
        m = T.ndim // 2; T = torch.diagonal(T, dim1=j, dim2=j + m).sum(-1)
    return T
def ratio(R):
    nn = torch.sum(torch.abs(R) ** 2); tot = 0
    for i in range(n):
        others = [j for j in range(n) if j != i]
        for r in range(len(others) + 1):
            for K in itertools.combinations(others, r):
                tot = tot + (-1.0 / d) ** r * torch.sum(torch.abs(ptr(R, [i] + list(K))) ** 2)
    return tot / d / nn
best = 0
for s in range(NS):
    torch.manual_seed(7 + s)
    A = torch.randn(D, 2, dtype=torch.complex128, requires_grad=True); B = torch.randn(D, 2, dtype=torch.complex128, requires_grad=True)
    o = torch.optim.Adam([A, B], lr=0.02)
    for it in range(2500):
        v = ratio(A @ torch.conj_physical(B).T); o.zero_grad(); (-v).backward(); o.step()
    best = max(best, v.item()); print(f"start {s}: ratio {v.item():.6f}", flush=True)
print(f"n={n} d={d}: max sum_i W_[n]-i / ||R||^2 = {best:.6f}   (FO1 bound 2/d = {2/d:.6f})")
