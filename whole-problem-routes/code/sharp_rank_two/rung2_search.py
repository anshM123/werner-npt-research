# Rung 2 of the Bernstein ladder:  R2(C) = S2 - 2(n-1) S1 + ((n-1)(3n+2)/2) S0 >= 0 for rank-2 C,
# S_k = sum_{|J|=k} ||Tr_J C||^2.  Adversarial minimisation of R2/S0 over complex rank-2 C on (C^3)^n.
import torch, itertools, sys, math
torch.set_default_dtype(torch.float64)
n = int(sys.argv[1]); runs = int(sys.argv[2]); d = 3
subsets = [J for J in itertools.product([0,1], repeat=n) if sum(J) <= 2]
def norms(Z):
    S = [0.0, 0.0, 0.0]
    for J in subsets:
        cur = Z.reshape([d]*(2*n)); m = n
        for i in sorted([i for i in range(n) if J[i]], reverse=True):
            cur = torch.diagonal(cur, dim1=i, dim2=i+m).sum(-1); m -= 1
        S[sum(J)] = S[sum(J)] + (cur.abs()**2).sum()
    return S
best = 1e9
for r in range(runs):
    torch.manual_seed(500 + 17*n + r)
    Ar = torch.randn(d**n, 2, requires_grad=True); Ai = torch.randn(d**n, 2, requires_grad=True)
    Br = torch.randn(d**n, 2, requires_grad=True); Bi = torch.randn(d**n, 2, requires_grad=True)
    opt = torch.optim.Adam([Ar, Ai, Br, Bi], lr=0.02)
    for it in range(3000):
        A = torch.complex(Ar, Ai); B = torch.complex(Br, Bi)
        Z = A @ B.conj().T
        S0, S1, S2 = norms(Z)
        val = (S2 - 2*(n-1)*S1 + (n-1)*(3*n+2)/2*S0)/S0
        opt.zero_grad(); val.backward(); opt.step()
    best = min(best, val.item())
    with torch.no_grad():
        sv = torch.linalg.svdvals(Z); sv = sv[:2]/sv[0]
    print(f"run {r}: R2/S0 = {val.item():+.3e}   singular value ratio {sv[1].item():.4f}", flush=True)
print("min R2/S0 =", best)
