"""Local minima of q_n(t, C)/||C||^2 over rank-two C at t < 1/2: value, local joint ranks, singular values, best cut."""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64)
n, runs, seed, t = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
d = 3 ** n
def q_t(C):
    T = C.reshape([3] * (2 * n)); X = T
    for s in range(n):
        tr = torch.diagonal(X, dim1=s, dim2=n + s).sum(-1)
        full = torch.einsum('ab,...->...ab', torch.eye(3, dtype=X.dtype), tr)
        X = X - t * torch.movedim(full, (2 * n - 2, 2 * n - 1), (s, n + s))
    return (T.conj() * X).sum().real
def joint(C, tol=1e-6):
    U, S, Vh = np.linalg.svd(C); out = []
    for s in range(n):
        blocks = [np.moveaxis(M.reshape([3]*n + [2]), s, 0).reshape(3, -1) for M in (U[:, :2], Vh[:2].conj().T)]
        Tj = np.concatenate(blocks, 1); ev = np.linalg.eigvalsh(Tj @ Tj.conj().T); out.append(int((ev > tol*ev.max()).sum()))
    return out, S[:2] / S[0]
sc = (1 - 2*t) * (1 - t) ** (n - 1)
rng = np.random.default_rng(seed)
for r in range(runs):
    ps = [torch.tensor(rng.standard_normal((d, 2)), requires_grad=True) for _ in range(4)]
    opt = torch.optim.LBFGS(ps, lr=1, max_iter=4000, tolerance_grad=1e-13, tolerance_change=1e-16, history_size=50, line_search_fn='strong_wolfe')
    def closure():
        opt.zero_grad(); A = torch.complex(ps[0], ps[1]); B = torch.complex(ps[2], ps[3]); C = A @ B.conj().T
        v = q_t(C) / (C.real**2 + C.imag**2).sum(); v.backward(); return v
    for _ in range(3): opt.step(closure)
    with torch.no_grad():
        A = torch.complex(ps[0], ps[1]); B = torch.complex(ps[2], ps[3]); C = A @ B.conj().T
        v = (q_t(C) / (C.real**2 + C.imag**2).sum()).item(); Cn = C.numpy()
    jr, sv = joint(Cn)
    print(f"run {r}: value/SC = {v/sc:.6f}  joint ranks {jr}  s2/s1 = {sv[1]:.4f}", flush=True)
