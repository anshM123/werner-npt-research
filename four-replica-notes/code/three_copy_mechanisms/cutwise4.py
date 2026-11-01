"""PSD rank-two endpoint at n = 4 via cut entropies (Wu-Zou style).  For orthonormal u, v in (C^3)^(x)4:
   g(x,y) = sum_T (-1)^{n-|T|} 2^{|T|} Tr(rho^x_T rho^y_T)   (= 2^n q-form; g(w,w) = 2^n q_4(P_w)),
   PSD endpoint  <=>  -g(u,v) <= sqrt(g(u,u) g(v,v)).
Candidate certificate (alpha + beta = 2):
   (I1)  -g(u,v) <= 1 + alpha sum_k sqrt(e_k(u) e_k(v)) + beta sum_c sqrt(e_c(u) e_c(v))   [k: 1|3 cuts, c: 2|2 cuts]
   (I2)  g(w,w)  >= 1 + alpha sum_k e_k(w) + beta sum_c e_c(w)
Then Cauchy-Schwarz gives the endpoint.  This script maximises the violation of I1 and of I2."""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64)
n, d = 4, 3
alpha = float(sys.argv[1]); beta = float(sys.argv[2])
which = sys.argv[3] if len(sys.argv) > 3 else 'I1'
starts = int(sys.argv[4]) if len(sys.argv) > 4 else 10
seed = int(sys.argv[5]) if len(sys.argv) > 5 else 0
subsets = [T for k in range(n + 1) for T in itertools.combinations(range(n), k)]
cuts1 = [(k,) for k in range(n)]
cuts2 = [(0, 1), (0, 2), (0, 3)]
def red(x, T):
    X = x.reshape((d,) * n)
    rest = [i for i in range(n) if i not in T]
    X = X.permute(*T, *rest).reshape(d ** len(T), d ** len(rest))
    return X @ X.conj().T
def tT(x, y, T):
    if len(T) == 0: return torch.real(torch.vdot(x, x) * torch.vdot(y, y))
    return torch.real(torch.trace(red(x, T) @ red(y, T)))
def g(x, y):
    return sum(((-1) ** (n - len(T))) * (2 ** len(T)) * tT(x, y, T) for T in subsets)
def ent(x, T):
    return 1 - tT(x, x, T)
def ortho(a, b):
    u = a / torch.linalg.norm(a); b = b - torch.vdot(u, b) * u
    return u, b / torch.linalg.norm(b)
rng = np.random.default_rng(seed); D = d ** n
best = -1e9
for s in range(starts):
    if which == 'I1':
        z = torch.tensor(rng.standard_normal((4, D)) * (1 if s % 2 else 0.3), requires_grad=True)
    else:
        z = torch.tensor(rng.standard_normal((2, D)), requires_grad=True)
    opt = torch.optim.Adam([z], lr=0.02)
    for it in range(2500):
        if which == 'I1':
            u, v = ortho(torch.complex(z[0], z[1]), torch.complex(z[2], z[3]))
            rhs = 1 + alpha * sum(torch.sqrt(torch.clamp(ent(u, T) * ent(v, T), min=1e-30)) for T in cuts1) \
                    + beta * sum(torch.sqrt(torch.clamp(ent(u, T) * ent(v, T), min=1e-30)) for T in cuts2)
            viol = -g(u, v) - rhs
        else:
            w = torch.complex(z[0], z[1]); w = w / torch.linalg.norm(w)
            viol = 1 + alpha * sum(ent(w, T) for T in cuts1) + beta * sum(ent(w, T) for T in cuts2) - g(w, w)
        opt.zero_grad(); (-viol).backward(); opt.step()
    best = max(best, viol.item())
    print(f"{which} alpha={alpha} beta={beta} start {s}: violation {viol.item():+.4e}", flush=True)
print(f"{which} alpha={alpha} beta={beta}: MAX violation {best:+.4e}", flush=True)
