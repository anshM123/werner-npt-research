"""Section 9, the n-site version (9.4') of (9.4) at rank two (a conjecture):
   sum_s ||Tr_s C||^2 <= (n-1)||C||^2 + (1/2)||C||_1^2 + (1/2)|Tr C|^2   for rank(C) <= 2.
Random-start minimisation (L-BFGS) of [RHS - LHS]/||C||^2 over complex factors C = A B^*.
usage: python multisite_refined.py runs seed d_1 ... d_n"""
import sys, numpy as np, torch
torch.set_default_dtype(torch.float64)
runs, seed = int(sys.argv[1]), int(sys.argv[2]); dims = [int(x) for x in sys.argv[3:]]; n = len(dims); D = int(np.prod(dims))
def gap(p):
    A = torch.complex(p[0], p[1]); B = torch.complex(p[2], p[3]); C = A @ B.conj().T
    n2 = (C.real ** 2 + C.imag ** 2).sum()
    GA, GB = A.conj().T @ A, B.conj().T @ B
    s1s2 = torch.sqrt(torch.clamp(((GA[0, 0] * GA[1, 1] - GA[0, 1] * GA[1, 0]) * (GB[0, 0] * GB[1, 1] - GB[0, 1] * GB[1, 0])).real, min=1e-24))
    tr = torch.diagonal(C).sum(); T = C.reshape(dims + dims); lhs = 0.
    for s in range(n):
        X = torch.diagonal(T, dim1=s, dim2=n + s).sum(-1); lhs = lhs + (X.real ** 2 + X.imag ** 2).sum()
    return ((n - 1) * n2 + 0.5 * (n2 + 2 * s1s2) + 0.5 * (tr.real ** 2 + tr.imag ** 2) - lhs) / n2, s1s2 / n2
rng = np.random.default_rng(seed); best, at = np.inf, None
for r in range(runs):
    p = [torch.tensor(rng.standard_normal((D, 2)), requires_grad=True) for _ in range(4)]
    opt = torch.optim.LBFGS(p, lr=1, max_iter=2000, tolerance_grad=1e-13, tolerance_change=1e-16, history_size=50, line_search_fn='strong_wolfe')
    def cl():
        opt.zero_grad(); v = gap(p)[0]; v.backward(); return v
    for _ in range(3): opt.step(cl)
    with torch.no_grad(): v, ratio = gap(p)
    if v.item() < best: best, at = v.item(), ratio.item()
print(f"dims={dims} runs={runs} seed={seed}: min [(n-1)|C|^2 + |C|_1^2/2 + |TrC|^2/2 - sum_s |Tr_s C|^2]/|C|^2 = {best:.3e}   (s1 s2/|C|^2 there: {at:.3f})")
