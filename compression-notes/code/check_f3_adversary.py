"""Section 11, F3: primal adversary.  The filters (A, B) : C^{3^M} -> C^3 take ascent steps on the distance of the
target X = (A (x) B) rho^{(x)M} (A (x) B)^* / Tr from K_1(3) = cone{products, (a (x) b) rho (a (x) b)^*}, estimated by a
warm-started fit with J filterings of rho plus P products (envelope gradient: the fit is held fixed during the step,
then refreshed).  A penalty keeps the target NPT (partial-transpose minimum at most -0.02); the worst NPT target is
refitted from fresh random starts at the end.  The run starts from the single-copy filtering (all copies but one
contracted) moved by delta0 toward random filters (delta0 < 0: random start).
usage: python check_f3_adversary.py M J P outer seed delta0   (e.g. 3 16 80 60 0 0.2; needs PyTorch, a few minutes)
"""
import sys, numpy as np, torch
torch.set_default_dtype(torch.float64)
M, J, P, OUTER, SEED = (int(a) for a in sys.argv[1:6]); DELTA0 = float(sys.argv[6]); m, D = 1, 3
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3): F[3 * i + j, 3 * j + i] = 1
rho = np.eye(9) - F / 2
def power(k):
    R = np.eye(1)
    for _ in range(k): R = np.kron(R, rho)
    perm = [2 * i for i in range(k)] + [2 * i + 1 for i in range(k)]
    return R.reshape([3] * (4 * k)).transpose(perm + [2 * k + q for q in perm]).reshape(9 ** k, 9 ** k)
Rm = torch.tensor(power(m), dtype=torch.complex128)
din, din1 = 3 ** m, 3 ** M
import itertools
subsets = [S for r in range(M + 1) for S in itertools.combinations(range(M), r)]
def filtM(A, B):
    At = A.reshape([D] + [3] * M); Bt = B.reshape([D] + [3] * M); Y = 0
    for Jn in subsets:
        Jc = [s_ for s_ in range(M) if s_ not in Jn]
        pa = [0] + [1 + s_ for s_ in Jn] + [1 + s_ for s_ in Jc]
        Ap = At.permute(pa).reshape(D, 3 ** len(Jn), 3 ** len(Jc)); Bp = Bt.permute(pa).reshape(D, 3 ** len(Jn), 3 ** len(Jc))
        M1 = torch.einsum('axc,eyc->aexy', Ap, Ap.conj()); M2 = torch.einsum('byc,fxc->bfyx', Bp, Bp.conj())
        Y = Y + (-0.5) ** len(Jn) * torch.einsum('aexy,bfyx->abef', M1, M2).reshape(D * D, D * D)
    return Y
torch.manual_seed(SEED)
def target(A, B):
    X = filtM(A, B); X = (X + X.conj().T) / 2
    return X / torch.trace(X).real
def model(p):
    V = torch.einsum('ka,kb->kab', p[2], p[3]).reshape(p[2].shape[0], D * D)
    Kt = torch.einsum('jab,jcd->jacbd', p[0], p[1]).reshape(p[0].shape[0], D * D, din * din)
    return V.T @ V.conj() + (Kt @ Rm @ Kt.conj().transpose(1, 2)).sum(0)
def refit(p, X, iters):
    opt = torch.optim.LBFGS(p, lr=1, max_iter=iters, history_size=50, line_search_fn='strong_wolfe', tolerance_grad=1e-15, tolerance_change=1e-18)
    def closure():
        opt.zero_grad(); loss = torch.linalg.norm(model(p) - X) ** 2; loss.backward(); return loss
    opt.step(closure)
    with torch.no_grad(): return float(torch.linalg.norm(model(p) - X) / torch.linalg.norm(X))
def pt(Y): return Y.reshape(D, D, D, D).permute(0, 3, 2, 1).reshape(D * D, D * D)
if DELTA0 >= 0:
    u = torch.zeros(3 ** (M - 1), dtype=torch.complex128); u[0] = 1
    A0 = torch.einsum('ax,y->axy', torch.eye(3, dtype=torch.complex128), u).reshape(3, din1)
    A = (A0 / torch.linalg.norm(A0) + DELTA0 * (lambda R: R / torch.linalg.norm(R))(torch.randn(D, din1, dtype=torch.complex128))).requires_grad_()
    B = (A0 / torch.linalg.norm(A0) + DELTA0 * (lambda R: R / torch.linalg.norm(R))(torch.randn(D, din1, dtype=torch.complex128))).requires_grad_()
else:
    A = torch.randn(D, din1, dtype=torch.complex128).requires_grad_(); B = torch.randn(D, din1, dtype=torch.complex128).requires_grad_()
p = [(0.4 * torch.randn(*s, dtype=torch.complex128)).requires_grad_() for s in ((J, D, din), (J, D, din), (P, D), (P, D))]
with torch.no_grad(): X = target(A, B)
r = refit(p, X, 3000)
optA = torch.optim.Adam([A, B], lr=0.01)
best = (0.0, X.detach().clone())
for it in range(OUTER):
    X = target(A, B)
    with torch.no_grad(): Mp = model(p)
    res = torch.linalg.norm(Mp - X) / torch.linalg.norm(X)
    ptm_t = torch.linalg.eigvalsh(pt(X))[0]
    obj = res - 10.0 * torch.relu(ptm_t + 0.02)          # stay NPT (PT-min <= -0.02)
    optA.zero_grad(); (-obj).backward(); optA.step()
    with torch.no_grad(): X = target(A, B)
    r = refit(p, X.detach(), 300)
    if r > best[0] and float(torch.linalg.eigvalsh(pt(X.detach()))[0]) < -0.005: best = (r, X.detach().clone())
    if it % 10 == 0 or it == OUTER - 1:
        ptm = float(torch.linalg.eigvalsh(pt(X.detach()))[0])
        print(f"F3 M = {M}, seed {SEED}, step {it}: residual {r:.3e} (PT-min {ptm:+.4f}); best so far {best[0]:.3e}", flush=True)
# re-check the worst target from fresh random starts, with a larger model
Xw = best[1]
def fresh(X, J, P, RS=3):
    out = np.inf
    for _ in range(RS):
        q = [(0.4 * torch.randn(*s, dtype=torch.complex128)).requires_grad_() for s in ((J, D, din), (J, D, din), (P, D), (P, D))]
        out = min(out, refit(q, X, 4000))
    return out
print(f"F3 M = {M}, seed {SEED}: worst target re-fitted from fresh starts: {J} + {P}: {fresh(Xw, J, P):.3e}; "
      f"{2 * J} + {2 * P}: {fresh(Xw, 2 * J, 2 * P):.3e}", flush=True)
np.save(f"adv_f3_M{M}_s{SEED}.npy", Xw.numpy())
