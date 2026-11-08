"""Section 9, Conjecture FB (merging two sites): fusion budget for the trace-norm gap.  Delta_dims(C) = q_dims(t,C) - (1-t)^(N-1)(||C||^2 - t||C||_1^2).
Test  Delta_(d1,d2,rest)(C) >= gamma * Delta_(d1 d2, rest)(C)  for rank-two C: minimise  Delta_full - gamma*Delta_merged.
usage: python fusion_budget.py t gamma runs seed d1 d2 ..."""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64)
t, gam, runs, seed = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
dims = [int(x) for x in sys.argv[5:]]
mdims = [dims[0] * dims[1]] + dims[2:]
D = int(np.prod(dims))
def q(C, dd):
    n = len(dd); T = C.reshape(dd + dd); tot = 0.
    for r in range(n + 1):
        for K in itertools.combinations(range(n), r):
            X, rem = T, list(range(n))
            for s in sorted(K, reverse=True):
                i = rem.index(s); X = torch.diagonal(X, dim1=i, dim2=len(rem) + i).sum(-1); rem.pop(i)
            tot = tot + (-t) ** len(K) * (X.real ** 2 + X.imag ** 2).sum()
    return tot
def parts(p):
    A = torch.complex(p[0], p[1]); B = torch.complex(p[2], p[3]); C = A @ B.conj().T
    n2 = (C.real ** 2 + C.imag ** 2).sum()
    GA, GB = A.conj().T @ A, B.conj().T @ B
    s1s2 = torch.sqrt(torch.clamp(((GA[0,0]*GA[1,1]-GA[0,1]*GA[1,0])*(GB[0,0]*GB[1,1]-GB[0,1]*GB[1,0])).real, min=1e-24))
    L = n2 - t * (n2 + 2 * s1s2)
    full = q(C, dims) - (1 - t) ** (len(dims) - 1) * L
    merged = q(C, mdims) - (1 - t) ** (len(mdims) - 1) * L
    return full, merged, n2
def obj(p):
    f, m, n2 = parts(p); return (f - gam * m) / n2
rng = np.random.default_rng(seed); best = 1e9; info = None
for r in range(runs):
    p = [torch.tensor(rng.standard_normal((D, 2)), requires_grad=True) for _ in range(4)]
    opt = torch.optim.LBFGS(p, lr=1, max_iter=2000, tolerance_grad=1e-13, tolerance_change=1e-16, history_size=50, line_search_fn='strong_wolfe')
    def cl():
        opt.zero_grad(); v = obj(p); v.backward(); return v
    for _ in range(3): opt.step(cl)
    with torch.no_grad():
        f, m, n2 = parts(p)
    v = ((f - gam * m) / n2).item()
    if v < best: best = v; info = (f.item() / n2.item(), m.item() / n2.item())
print(f"t={t} gamma={gam} dims={dims}->{mdims}: min (Delta_full - gamma Delta_merged)/|C|^2 = {best:.3e}   (Delta_full, Delta_merged) there = ({info[0]:.3e}, {info[1]:.3e})")
