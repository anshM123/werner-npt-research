"""Werner-line minimum: min over rank-two C of q_n(t, C)/||C||^2 versus the thin product value (1-2t)(1-t)^(n-1)."""
import sys, numpy as np, torch
torch.set_default_dtype(torch.float64)
n, runs, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
ts = [float(x) for x in sys.argv[4].split(',')]
d = 3 ** n
def q_t(C, t):
    T = C.reshape([3] * (2 * n)); X = T
    for s in range(n):
        tr = torch.diagonal(X, dim1=s, dim2=n + s).sum(-1)
        full = torch.einsum('ab,...->...ab', torch.eye(3, dtype=X.dtype), tr)
        X = X - t * torch.movedim(full, (2 * n - 2, 2 * n - 1), (s, n + s))
    return (T.conj() * X).sum().real
rng = np.random.default_rng(seed)
for t in ts:
    best = None
    for r in range(runs):
        ps = [torch.tensor(rng.standard_normal((d, 2)), requires_grad=True) for _ in range(4)]
        opt = torch.optim.LBFGS(ps, lr=1, max_iter=3000, tolerance_grad=1e-12, tolerance_change=1e-15, history_size=50, line_search_fn='strong_wolfe')
        def closure():
            opt.zero_grad()
            A = torch.complex(ps[0], ps[1]); B = torch.complex(ps[2], ps[3]); C = A @ B.conj().T
            v = q_t(C, t) / (C.real ** 2 + C.imag ** 2).sum(); v.backward(); return v
        for _ in range(2): opt.step(closure)
        with torch.no_grad():
            A = torch.complex(ps[0], ps[1]); B = torch.complex(ps[2], ps[3]); C = A @ B.conj().T
            v = (q_t(C, t) / (C.real ** 2 + C.imag ** 2).sum()).item()
        best = v if best is None else min(best, v)
    sc = (1 - 2 * t) * (1 - t) ** (n - 1)
    print(f"n = {n}, t = {t}: min found {best:.6f}   thin product (1-2t)(1-t)^(n-1) = {sc:.6f}   ratio {best/sc:.4f}", flush=True)
