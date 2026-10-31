"""Section 11: are PPT 3 x 3 local filterings of rho^{(x)m} separable?  F3 needs this for PPT-entangled states of
Tiles type (their range contains no product vector, so they lie outside every K_j).  The filters take ascent steps on
the residual of a warm-started fit by P product terms, with a penalty that keeps the target strictly PPT (PT-min at
least 0.002); the worst PPT target is refitted from fresh starts.  Control: the Tiles state with 2% white noise.
usage: python check_ppt_filterings.py m P outer seed   (e.g. 2 40 80 0; needs PyTorch, a few minutes)
"""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64)
m, P, OUTER, SEED = (int(a) for a in sys.argv[1:5]); D = 3
subsets = [S for r in range(m + 1) for S in itertools.combinations(range(m), r)]
def filt(A, B):
    At, Bt = A.reshape([D] + [3] * m), B.reshape([D] + [3] * m); Y = 0
    for S in subsets:
        Sc = [s for s in range(m) if s not in S]
        pa = [0] + [1 + s for s in S] + [1 + s for s in Sc]
        Ap = At.permute(pa).reshape(D, 3 ** len(S), 3 ** len(Sc)); Bp = Bt.permute(pa).reshape(D, 3 ** len(S), 3 ** len(Sc))
        M1 = torch.einsum('axc,eyc->aexy', Ap, Ap.conj()); M2 = torch.einsum('byc,fxc->bfyx', Bp, Bp.conj())
        Y = Y + (-0.5) ** len(S) * torch.einsum('aexy,bfyx->abef', M1, M2).reshape(D * D, D * D)
    Y = (Y + Y.conj().T) / 2
    return Y / torch.trace(Y).real
def pt(Y): return Y.reshape(3, 3, 3, 3).permute(0, 3, 2, 1).reshape(9, 9)
def model(p):
    V = torch.einsum('ka,kb->kab', p[0], p[1]).reshape(p[0].shape[0], 9)
    return V.T @ V.conj()
def refit(p, X, iters):
    opt = torch.optim.LBFGS(p, lr=1, max_iter=iters, history_size=50, line_search_fn='strong_wolfe', tolerance_grad=1e-15, tolerance_change=1e-18)
    def closure():
        opt.zero_grad(); loss = torch.linalg.norm(model(p) - X) ** 2; loss.backward(); return loss
    opt.step(closure)
    with torch.no_grad(): return float(torch.linalg.norm(model(p) - X) / torch.linalg.norm(X))
def fresh(X, P, RS=3):
    out = np.inf
    for _ in range(RS):
        q = [(0.3 * torch.randn(P, 3, dtype=torch.complex128)).requires_grad_() for _ in range(2)]
        out = min(out, refit(q, X, 4000))
    return out
# control: Tiles UPB state
v = [np.array(x, dtype=complex) for x in ([1, 0, 0], [0, 1, 0], [0, 0, 1])]
e0, e1, e2 = v
upb = [np.kron(e0, (e0 - e1) / 2 ** .5), np.kron(e2, (e1 - e2) / 2 ** .5), np.kron((e0 - e1) / 2 ** .5, e2),
       np.kron((e1 - e2) / 2 ** .5, e0), np.kron((e0 + e1 + e2) / 3 ** .5, (e0 + e1 + e2) / 3 ** .5)]
Pu = sum(np.outer(x, x.conj()) for x in upb)
T = (np.eye(9) - Pu) / 4; T = 0.98 * T + 0.02 * np.eye(9) / 9
torch.manual_seed(SEED)
print(f"control Tiles + 2% noise: PT-min {np.linalg.eigvalsh(pt(torch.tensor(T)).numpy())[0]:+.4f}, product fit ({P} terms) "
      f"{fresh(torch.tensor(T, dtype=torch.complex128), P):.2e}", flush=True)
A = torch.randn(D, 3 ** m, dtype=torch.complex128).requires_grad_(); B = torch.randn(D, 3 ** m, dtype=torch.complex128).requires_grad_()
p = [(0.3 * torch.randn(P, 3, dtype=torch.complex128)).requires_grad_() for _ in range(2)]
with torch.no_grad(): X = filt(A, B)
refit(p, X, 3000)
optA = torch.optim.Adam([A, B], lr=0.01)
best = (0.0, X.detach().clone())
for it in range(OUTER):
    X = filt(A, B)
    with torch.no_grad(): Mp = model(p)
    res = torch.linalg.norm(Mp - X) / torch.linalg.norm(X)
    ptm = torch.linalg.eigvalsh(pt(X))[0]
    obj = res - 10.0 * torch.relu(0.002 - ptm)
    optA.zero_grad(); (-obj).backward(); optA.step()
    with torch.no_grad(): X = filt(A, B)
    r = refit(p, X.detach(), 300)
    ptv = float(torch.linalg.eigvalsh(pt(X.detach()))[0])
    if r > best[0] and ptv > 0: best = (r, X.detach().clone())
    if it % 10 == 0 or it == OUTER - 1:
        print(f"PS m = {m}, seed {SEED}, step {it}: residual {r:.3e} (PT-min {ptv:+.4f}); best so far {best[0]:.3e}", flush=True)
print(f"PS m = {m}, seed {SEED}: worst PPT target refit from fresh starts: {P}: {fresh(best[1], P):.3e}; {2 * P}: "
      f"{fresh(best[1], 2 * P):.3e}", flush=True)
