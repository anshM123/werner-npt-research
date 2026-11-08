# Adversarial test of polydisk zero-freeness: F_Z(t) = sum_J t^J ||Tr_J Z||^2, Z rank 2 on (C^3)^n.
# F = a + b t_1 (affine in t_1).  Zero-free on |t_s|<1/2 for all s  <=>  |a|/|b| >= 1/2 whenever |t_s|<=1/2 (s>=2).
import torch, itertools, sys, math
torch.set_default_dtype(torch.float64)
n = int(sys.argv[1]); runs = int(sys.argv[2]); mode = sys.argv[3] if len(sys.argv) > 3 else 'multi'
d = 3
subsets = list(itertools.product([0,1], repeat=n))
def trace_norms(Z):
    out = {}
    for J in subsets:
        cur = Z.reshape([d]*(2*n)); m = n
        for i in sorted([i for i in range(n) if J[i]], reverse=True):
            cur = torch.diagonal(cur, dim1=i, dim2=i+m).sum(-1); m -= 1
        out[J] = (cur.abs()**2).sum()
    return out
best = 1e9
for r in range(runs):
    torch.manual_seed(1000*n + r)
    Ar = torch.randn(d**n, 2, requires_grad=True); Ai = torch.randn(d**n, 2, requires_grad=True)
    Br = torch.randn(d**n, 2, requires_grad=True); Bi = torch.randn(d**n, 2, requires_grad=True)
    rad = torch.randn(n, requires_grad=True); ang = (2*math.pi*torch.rand(n)).requires_grad_()
    opt = torch.optim.Adam([Ar, Ai, Br, Bi, rad, ang], lr=0.03)
    for it in range(2500):
        A = torch.complex(Ar, Ai); B = torch.complex(Br, Bi)
        Z = A @ B.conj().T
        c = trace_norms(Z)
        if mode == 'multi':
            t = 0.5*torch.sigmoid(rad)*torch.exp(1j*ang)
        else:
            t = (0.5*torch.sigmoid(rad[0])*torch.exp(1j*ang[0])).repeat(n)
        a = 0; b = 0
        for J, v in c.items():
            w = 1
            for s in range(1, n):
                if J[s]: w = w * t[s]
            if J[0]: b = b + w*v
            else: a = a + w*v
        if mode == 'multi':
            ratio = a.abs()/b.abs()      # need >= 1/2 ; t_1 is free
        else:
            # univariate: F(t) = sum_k t^k m_k ; check |F(t)| small relative to scale at t on disk
            F = a + t[0]*b
            ratio = F.abs()/ (sum(c.values()))
        opt.zero_grad(); ratio.backward(); opt.step()
    best = min(best, ratio.item())
    print(f"run {r}: ratio {ratio.item():.6f}  |t| = {[round(x,3) for x in (0.5*torch.sigmoid(rad)).tolist()]}", flush=True)
print("min ratio", best, "(multi: need >= 0.5)")
