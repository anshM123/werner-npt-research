# Polydisk test restricted to BALANCED rank-2 projectors (real), excluding product zeros.
# P = Q Q^T, Q real 3^n x 2 orthonormal; penalty forces each one-site marginal of P/2 to be I/3.
# Objective: min over P and t' (|t_s|<=1/2, s>=2) of |a|/|b|, F = a + b t_1.
import torch, itertools, sys, math
torch.set_default_dtype(torch.float64)
n = int(sys.argv[1]); runs = int(sys.argv[2]); lam = float(sys.argv[3]) if len(sys.argv)>3 else 50.0
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
def marg_pen(P):
    pen = 0
    T = P.reshape([d]*(2*n))
    for s in range(n):
        cur = T; m = n
        for i in sorted([i for i in range(n) if i != s], reverse=True):
            cur = torch.diagonal(cur, dim1=i, dim2=i+m).sum(-1); m -= 1
        pen = pen + ((cur/2 - torch.eye(d)/3)**2).sum()
    return pen
best = 1e9
for r in range(runs):
    torch.manual_seed(77*n + r)
    X = torch.randn(d**n, 2, requires_grad=True)
    rad = torch.randn(n, requires_grad=True); ang = (2*math.pi*torch.rand(n)).requires_grad_()
    opt = torch.optim.Adam([X, rad, ang], lr=0.02)
    for it in range(3000):
        Q, _ = torch.linalg.qr(X); P = Q @ Q.T
        c = trace_norms(P.to(torch.complex128))
        t = 0.5*torch.sigmoid(rad)*torch.exp(1j*ang)
        a = 0; b = 0
        for J, v in c.items():
            w = 1
            for s in range(1, n):
                if J[s]: w = w * t[s]
            if J[0]: b = b + w*v
            else: a = a + w*v
        ratio = a.abs()/b.abs()
        pen = marg_pen(P)
        loss = ratio + lam*pen
        opt.zero_grad(); loss.backward(); opt.step()
    best = min(best, ratio.item())
    print(f"run {r}: ratio {ratio.item():.6f} pen {pen.item():.2e} |t| = {[round(x,3) for x in (0.5*torch.sigmoid(rad)).tolist()]} arg = {[round(x % (2*math.pi),2) for x in ang.tolist()]}", flush=True)
print("min ratio", best)
