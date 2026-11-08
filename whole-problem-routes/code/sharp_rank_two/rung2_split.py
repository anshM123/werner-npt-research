# Rung 2 split.  e_a = ||Tr_a C||^2 - ||C||^2,  f_ab = ||Tr_ab C||^2 - ||Tr_a C||^2 - ||Tr_b C||^2 + ||C||^2.
# (P1): sum_{a in A} e_a <= ||C||^2 for every subset A (=> sum_a e_a^+ <= ||C||^2).
# (P2): f_ab + e_a^- + e_b^- >= 0 (other copies untraced).   Both for rank-2 C on (C^3)^n.
import torch, itertools, sys, math
torch.set_default_dtype(torch.float64)
mode = sys.argv[1]; n = int(sys.argv[2]); runs = int(sys.argv[3]); d = 3
def pt_norm2(Z, J):
    cur = Z.reshape([d]*(2*n)); m = n
    for i in sorted(J, reverse=True):
        cur = torch.diagonal(cur, dim1=i, dim2=i+m).sum(-1); m -= 1
    return (cur.abs()**2).sum()
def make(Ar, Ai, Br, Bi):
    return torch.complex(Ar, Ai) @ torch.complex(Br, Bi).conj().T
results = []
targets = ([tuple(A) for k in range(2, n+1) for A in itertools.combinations(range(n), k)] if mode == 'P1' else [(0, 1)])
for A in targets:
    best = -1e9 if mode == 'P1' else 1e9
    for r in range(runs):
        torch.manual_seed(7000 + 31*r + 7*len(A))
        P = [torch.randn(d**n, 2, requires_grad=True) for _ in range(4)]
        opt = torch.optim.Adam(P, lr=0.02)
        for it in range(2500):
            Z = make(*P); S0 = (Z.abs()**2).sum()
            if mode == 'P1':
                val = sum(pt_norm2(Z, [a]) - S0 for a in A)/S0
                loss = -val
            else:
                a, b = A
                ea = (pt_norm2(Z, [a]) - S0)/S0; eb = (pt_norm2(Z, [b]) - S0)/S0
                fab = (pt_norm2(Z, [a, b]) - pt_norm2(Z, [a]) - pt_norm2(Z, [b]) + S0)/S0
                beta = 60.0
                val = fab + torch.nn.functional.softplus(-ea, beta=beta) + torch.nn.functional.softplus(-eb, beta=beta)
                loss = val
            opt.zero_grad(); loss.backward(); opt.step()
        with torch.no_grad():
            if mode == 'P2':
                exact = (fab + torch.clamp(-ea, min=0) + torch.clamp(-eb, min=0)).item()
                best = min(best, exact)
            else:
                best = max(best, val.item())
    results.append((A, best))
    print(mode, 'n =', n, 'set', A, ('max sum e_a/||C||^2' if mode == 'P1' else 'min f_ab + e_a^- + e_b^-  (/||C||^2)'), f'{best:+.6f}', flush=True)
