"""OS_n test with autograd: maximise g(u,v) = sum_T (-1)^{n-|T|} 2^{|T|} Tr(rho^u_T rho^v_T) over orthonormal u, v
(T = empty contributes (-1)^n, T = [n] contributes 0 since u perp v).  OS_n: g(u,v) <= 1."""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64)
n = int(sys.argv[1]); NS = int(sys.argv[2]); seed = int(sys.argv[3]); iters = int(sys.argv[4]) if len(sys.argv) > 4 else 3000
d = 3; D = d ** n
subsets = [T for k in range(1, n) for T in itertools.combinations(range(n), k)]
def red(x, T):
    X = x.reshape((d,) * n); rest = [i for i in range(n) if i not in T]
    X = X.permute(*T, *rest).reshape(d ** len(T), d ** len(rest)); return X @ X.conj().T
def g(u, v):
    out = torch.tensor(float((-1) ** n))
    for T in subsets:
        out = out + ((-1) ** (n - len(T))) * (2 ** len(T)) * torch.real(torch.trace(red(u, T) @ red(v, T)))
    return out
def ortho(a, b):
    u = a / torch.linalg.norm(a); b = b - torch.vdot(u, b) * u; return u, b / torch.linalg.norm(b)
rng = np.random.default_rng(seed); best = -1e9
for s in range(NS):
    z = torch.tensor(rng.standard_normal((4, D)) * (1.0 if s % 2 else 0.3), requires_grad=True)
    opt = torch.optim.Adam([z], lr=0.02)
    for it in range(iters):
        u, v = ortho(torch.complex(z[0], z[1]), torch.complex(z[2], z[3]))
        val = g(u, v)
        opt.zero_grad(); (-val).backward(); opt.step()
    best = max(best, val.item())
    print(f"n={n} start {s}: g(u,v) = {val.item():.8f}", flush=True)
print(f"n={n}: MAX g(u,v) over {NS} starts = {best:.8f}  (OS_n: <= 1)", flush=True)
