"""Section 4: Perron-operator test for rank-two MPS operators of bond dimension D.
Transfer channel  T(X) = sum_{u,v} M_uv X M_uv^* - (1/2) sum_{u,v} M_uv X M_vu^*,  M_uv = A^u (x) B^v  (rho = I - F/2).
X_* = Perron operator (power iteration from I, trace-normalised).  If some Schmidt-rank-two omega has
<omega| X_*^Gamma |omega> < 0, then q_n(1/2, C) < 0 for every large n (suitable boundary vectors).
eta_D = min over isometries U, V : C^2 -> C^D of lambda_min((U (x) V)^* X_*^Gamma (U (x) V)).
Joint minimisation of eta over (A, B, U, V); periodic careful inner minimisation over (U, V) for the current X_*.
usage: python check_perron_operator.py D seed steps [iters]   (needs PyTorch; D = 5: ~10 min, D = 9: ~40 min for 1500 steps)"""
import sys, numpy as np, torch
torch.set_default_dtype(torch.float64)
D, SEED, STEPS = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); IT = int(sys.argv[4]) if len(sys.argv) > 4 else 120
torch.manual_seed(SEED)
def perron(A, B, iters=IT):
    M = [[torch.kron(A[u], B[v]) for v in range(3)] for u in range(3)]
    X = torch.eye(D * D, dtype=torch.complex128) / (D * D)
    for k in range(iters):
        Y = 0
        for u in range(3):
            for v in range(3):
                Y = Y + M[u][v] @ X @ M[u][v].conj().T - 0.5 * M[u][v] @ X @ M[v][u].conj().T
        Y = (Y + Y.conj().T) / 2
        Xn = Y / torch.trace(Y).real
        if k == iters - 1: conv = torch.linalg.norm(Xn - X).item()
        X = Xn
    return X, conv
def pt(X): return X.reshape(D, D, D, D).permute(0, 3, 2, 1).reshape(D * D, D * D)
def iso(W):
    Q, R = torch.linalg.qr(W); return Q
def eta_of(XG, U, V):
    K = torch.kron(iso(U), iso(V))
    return torch.linalg.eigvalsh(K.conj().T @ XG @ K)[0]
def inner_eta(XG, restarts=12, iters=300):
    best = np.inf
    for r in range(restarts):
        U = torch.randn(D, 2, dtype=torch.complex128, requires_grad=True); V = torch.randn(D, 2, dtype=torch.complex128, requires_grad=True)
        o = torch.optim.Adam([U, V], lr=0.05)
        for it in range(iters):
            e = eta_of(XG, U, V); o.zero_grad(); e.backward(); o.step()
        best = min(best, e.item())
    return best
A = (torch.randn(3, D, D, dtype=torch.complex128) / np.sqrt(D)).requires_grad_()
B = (torch.randn(3, D, D, dtype=torch.complex128) / np.sqrt(D)).requires_grad_()
U = torch.randn(D, 2, dtype=torch.complex128, requires_grad=True); V = torch.randn(D, 2, dtype=torch.complex128, requires_grad=True)
opt = torch.optim.Adam([A, B, U, V], lr=0.01)
best = np.inf
for st in range(STEPS):
    X, conv = perron(A, B)
    XG = pt(X)
    e = eta_of(XG, U, V)
    opt.zero_grad(); e.backward(); opt.step()
    if st % 50 == 0 or st == STEPS - 1:
        with torch.no_grad():
            Xd, conv = perron(A.detach(), B.detach(), iters=3 * IT)
            ptmin = torch.linalg.eigvalsh(pt(Xd))[0].item()
        ie = inner_eta(pt(Xd).detach())
        best = min(best, ie)
        print(f"D = {D}, seed {SEED}, step {st}: joint eta {e.item():+.3e}; careful eta(X_*) {ie:+.3e} (PT-min of X_* {ptmin:+.3e}, "
              f"power-iteration change {conv:.1e}); best {best:+.3e}", flush=True)
        if ie < -1e-6:
            np.savez(f"perron_neg_D{D}_s{SEED}_st{st}.npz", A=A.detach().numpy(), B=B.detach().numpy(), eta=ie)
np.savez(f"perron_final_D{D}_s{SEED}.npz", A=A.detach().numpy(), B=B.detach().numpy(), best=best)
