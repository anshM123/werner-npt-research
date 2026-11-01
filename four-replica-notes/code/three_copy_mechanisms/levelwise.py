"""Level-wise generalisation of the Wu-Zou exterior-square inequality (arXiv:2608.02647, Thm 1.2).
For orthonormal u, v in (C^d)^(x)n put p_K(x(x)y) = ||Pi_K (x (x) y)||^2, Pi_K = (x)_{i in K} P^-_i (x)_{i not in K} P^+_i,
p_k = sum_{|K|=k} p_K(u (x) v), A_k(w) = sum_{|K|=k} p_K(w (x) w).
Claim L_j:  3 p_{2j+1}(u,v) - p_{2j}(u,v) <= sqrt(A_{2j}(u) A_{2j}(v))   (n = 3, j = 1 is Wu-Zou Thm 1.2).
All L_j  =>  PSD rank-two endpoint at n copies (Cauchy-Schwarz).
Walsh: p_K(x,y) = 2^-n sum_T (-1)^{|K cap T|} Tr(rho^x_T rho^y_T)."""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64)
n = int(sys.argv[1]); d = int(sys.argv[2]) if len(sys.argv) > 2 else 3
starts = int(sys.argv[3]) if len(sys.argv) > 3 else 20
seed = int(sys.argv[4]) if len(sys.argv) > 4 else 0
subsets = [T for k in range(n + 1) for T in itertools.combinations(range(n), k)]
def red(x, T):
    # reduced density matrix of the pure (unnormalised) vector x on sites T
    X = x.reshape((d,) * n)
    rest = [i for i in range(n) if i not in T]
    X = X.permute(*T, *rest).reshape(d ** len(T), d ** len(rest))
    return X @ X.conj().T
def tvec(x, y):
    out = []
    for T in subsets:
        if len(T) == 0: out.append(torch.real(torch.vdot(x, x) * torch.vdot(y, y))); continue
        out.append(torch.real(torch.trace(red(x, T) @ red(y, T))))
    return torch.stack(out)
# Walsh matrix: p_K = 2^-n sum_T (-1)^{|K cap T|} t_T, grouped by |K|
Wk = np.zeros((n + 1, len(subsets)))
for K in subsets:
    for j, T in enumerate(subsets):
        Wk[len(K), j] += (-1) ** len(set(K) & set(T)) / 2 ** n
Wk = torch.tensor(Wk)
def levels(u, v):
    pk = Wk @ tvec(u, v); Au = Wk @ tvec(u, u); Av = Wk @ tvec(v, v)
    return pk, Au, Av
def ortho(a, b):
    u = a / torch.linalg.norm(a)
    b = b - torch.vdot(u, b) * u
    return u, b / torch.linalg.norm(b)
rng = np.random.default_rng(seed)
D = d ** n
for j in range(1, (n - 1) // 2 + 1):
    best = -1e9
    for s in range(starts):
        z = torch.tensor(rng.standard_normal((4, D)), requires_grad=True)
        opt = torch.optim.Adam([z], lr=0.03)
        for it in range(1500):
            a = torch.complex(z[0], z[1]); b = torch.complex(z[2], z[3])
            u, v = ortho(a, b)
            pk, Au, Av = levels(u, v)
            viol = 3 * pk[2 * j + 1] - pk[2 * j] - torch.sqrt(torch.clamp(Au[2 * j] * Av[2 * j], min=1e-30))
            opt.zero_grad(); (-viol).backward(); opt.step()
        best = max(best, viol.item())
        print(f"n={n} d={d} j={j} start {s}: max violation L_j = {viol.item():+.3e}  (p={pk.detach().numpy().round(4)}, A_u={Au[2*j].item():.4f}, A_v={Av[2*j].item():.4f})", flush=True)
    print(f"n={n} d={d} j={j}: MAX over starts of [3p_{2*j+1} - p_{2*j} - sqrt(A A')] = {best:+.3e}", flush=True)
