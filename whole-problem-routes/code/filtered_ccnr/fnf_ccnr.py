# Independent check of the agent's conjecture S1 at n copies: for local compressions Y = (A (x) B)^* W^{(x)n} (A (x) B),
# A: C^r -> H, B: C^r' -> H, bring Y to filter normal form (both marginals ~ identity, Sinkhorn-type iteration) and
# compute ||R(Y~)||_1 / Tr Y~ (realignment / CCNR).  S1 claims <= 1.  Maximise adversarially.
import torch, sys, math
torch.set_default_dtype(torch.float64)
n = int(sys.argv[1]); r = int(sys.argv[2]); rp = int(sys.argv[3]); runs = int(sys.argv[4]); t = float(sys.argv[5]) if len(sys.argv) > 5 else 0.5
d = 3; D = d**n
# W^{(x)n} on (A side: n copies) (x) (B side: n copies), ordering [A_1..A_n, B_1..B_n]
I9 = torch.eye(d*d, dtype=torch.complex128)
F = torch.zeros(d*d, d*d, dtype=torch.complex128)
for i in range(d):
    for j in range(d):
        F[i*d+j, j*d+i] = 1
W1 = I9 - t*F                                  # on (A_s, B_s)
Wn = W1
for _ in range(n-1):
    Wn = torch.kron(Wn, W1)                    # ordering (A1 B1 A2 B2 ...)
# reorder to (A1..An, B1..Bn)
perm = [2*s for s in range(n)] + [2*s+1 for s in range(n)]
Wt = Wn.reshape([d]*(4*n)).permute(perm + [2*n + p for p in perm]).reshape(D*D, D*D)
def msqrt_inv(M):
    # any G with G^* M G = I : inverse-adjoint Cholesky factor (gauge-free gradient)
    M = (M + M.conj().T)/2
    L = torch.linalg.cholesky(M + 1e-13*torch.eye(M.shape[0], dtype=M.dtype))
    return torch.linalg.inv(L).conj().T
def ratio(A, B, iters=60):
    K = torch.kron(A, B)
    Y = K.conj().T @ Wt @ K                      # (r r') x (r r')
    for _ in range(iters):                       # filter normal form
        Ya = Y.reshape(r, rp, r, rp).diagonal(dim1=1, dim2=3).sum(-1)   # Tr_B -> r x r
        Ga = msqrt_inv(Ya)
        Gk = torch.kron(Ga, torch.eye(rp, dtype=Y.dtype)); Y = Gk.conj().T @ Y @ Gk
        Yb = Y.reshape(r, rp, r, rp).diagonal(dim1=0, dim2=2).sum(-1)   # Tr_A -> r' x r'
        Gb = msqrt_inv(Yb)
        Gk = torch.kron(torch.eye(r, dtype=Y.dtype), Gb); Y = Gk.conj().T @ Y @ Gk
    # realignment: Y_{(a b),(a' b')} -> R_{(a a'),(b b')}
    R = Y.reshape(r, rp, r, rp).permute(0, 2, 1, 3).reshape(r*r, rp*rp)
    return torch.linalg.svdvals(R).sum() / torch.trace(Y).real
best = 0
for run in range(runs):
    torch.manual_seed(1000 + run)
    Ar = torch.randn(D, r, requires_grad=True); Ai = torch.randn(D, r, requires_grad=True)
    Br = torch.randn(D, rp, requires_grad=True); Bi = torch.randn(D, rp, requires_grad=True)
    opt = torch.optim.Adam([Ar, Ai, Br, Bi], lr=0.02)
    for it in range(400):
        A = torch.complex(Ar, Ai); B = torch.complex(Br, Bi)
        val = ratio(A, B)
        opt.zero_grad(); (-val).backward(); opt.step()
    best = max(best, val.item())
    print(f"run {run}: filtered CCNR ratio {val.item():.6f}", flush=True)
print(f"n={n} r={r} r'={rp} t={t}: max ratio {best:.6f}")
