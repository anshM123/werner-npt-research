"""Section 9: positive semidefinite C of rank at most r on (C^3)^n; minimises [2^n q_n(C) - (2||C||^2 - (Tr C)^2)]/(Tr C)^2.
At n = 3 the minimum is 0 for r = 2 and negative from r = 3 on, so the rank enters even for positive operators.
usage: python psd_rank.py n r runs"""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64)
n = int(sys.argv[1]); rk = int(sys.argv[2]); runs = int(sys.argv[3])
d = 3 ** n; subsets = [K for r in range(n + 1) for K in itertools.combinations(range(n), r)]
def qn(C):
    T = C.reshape([3] * n + [3] * n); tot = 0.
    for K in subsets:
        X, rem = T, list(range(n))
        for s in sorted(K, reverse=True):
            i = rem.index(s); X = torch.diagonal(X, dim1=i, dim2=len(rem) + i).sum(-1); rem.pop(i)
        tot = tot + (-0.5) ** len(K) * (X.real ** 2 + X.imag ** 2).sum()
    return tot
def gap(p):
    G = torch.complex(p[0], p[1]); C = G @ G.conj().T
    n2 = (C.real ** 2 + C.imag ** 2).sum(); tr = torch.trace(C).real
    return (2 ** n * qn(C) - (2 * n2 - tr ** 2)) / tr ** 2
rng = np.random.default_rng(3); out = []
for r in range(runs):
    p = [torch.tensor(rng.standard_normal((d, rk)), requires_grad=True) for _ in range(2)]
    opt = torch.optim.LBFGS(p, lr=1, max_iter=2000, tolerance_grad=1e-12, tolerance_change=1e-15, history_size=50, line_search_fn='strong_wolfe')
    def cl():
        opt.zero_grad(); v = gap(p); v.backward(); return v
    for _ in range(2): opt.step(cl)
    out.append(gap(p).item())
print(f"n={n} PSD rank<={rk}: min [2^n q_n - (2||C||^2 - (TrC)^2)]/(TrC)^2 = {min(out):.4e}")
