"""Primal adversary for L_m(6): maximise over filters (A, B): C^{3^(m+1)} -> C^6 the distance of
   X = (A(x)B) rho^{(x)(m+1)} (A(x)B)^* / Tr from K_m(6) = cone{products, (a(x)b) rho^{(x)m} (a(x)b)^*}, estimated by a
   warm-started fit with J filterings + P products.  Envelope gradient: the fit parameters are held fixed while (A, B)
   take an ascent step on the residual, then the fit is refreshed.  m = 1 (L_1(6), which looks false) is the control.
   usage: python check_l2six_adversary_npt.py J P outer delta  (m = 2, started from the three-copy target of Section 9)"""
import sys, numpy as np, torch
torch.set_default_dtype(torch.float64)
J, P, OUTER = (int(a) for a in sys.argv[1:4]); DELTA = float(sys.argv[4]); m = 2; SEED = int(100 * DELTA)
D = 6
F = np.zeros((9, 9))
for i in range(3):
    for j in range(3): F[3 * i + j, 3 * j + i] = 1
rho = np.eye(9) - F / 2
def power(k):
    R = np.eye(1)
    for _ in range(k): R = np.kron(R, rho)
    perm = [2 * i for i in range(k)] + [2 * i + 1 for i in range(k)]
    return R.reshape([3] * (4 * k)).transpose(perm + [2 * k + q for q in perm]).reshape(9 ** k, 9 ** k)
Rm = torch.tensor(power(m), dtype=torch.complex128); Rm1 = torch.tensor(power(m + 1), dtype=torch.complex128)
din, din1 = 3 ** m, 3 ** (m + 1)
torch.manual_seed(SEED)
def target(A, B):
    K = torch.kron(A, B); X = K @ Rm1 @ K.conj().T; X = (X + X.conj().T) / 2
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
R2_ = torch.tensor(power(2), dtype=torch.complex128)
def ptt(Y): return Y.reshape(D, D, D, D).permute(0, 3, 2, 1).reshape(D * D, D * D)
torch.manual_seed(300)
A2 = torch.randn(D, 9, dtype=torch.complex128); B2 = torch.randn(D, 9, dtype=torch.complex128)
A2 = A2.clone().requires_grad_(); B2 = B2.clone().requires_grad_(); o2 = torch.optim.Adam([A2, B2], lr=0.005)
for thr in (-0.01, -0.03):
    for it in range(4000):
        K = torch.kron(A2, B2); Yt = K @ R2_ @ K.conj().T; Yt = (Yt + Yt.conj().T) / 2; Yt = Yt / torch.trace(Yt).real
        f = torch.linalg.eigvalsh(ptt(Yt))[0]
        if f.item() < thr: break
        o2.zero_grad(); f.backward(); o2.step()
A2, B2 = A2.detach(), B2.detach()
g = torch.Generator().manual_seed(1)
u = torch.randn(3, dtype=torch.complex128, generator=g); u = u / torch.linalg.norm(u)
v = torch.randn(3, dtype=torch.complex128, generator=g); v = v / torch.linalg.norm(v)
RA = torch.randn(D, 27, dtype=torch.complex128, generator=g); RB = torch.randn(D, 27, dtype=torch.complex128, generator=g)
A = (torch.kron(A2, u.conj()[None, :]) + DELTA * RA * torch.linalg.norm(A2) / torch.linalg.norm(RA)).requires_grad_()
B = (torch.kron(B2, v.conj()[None, :]) + DELTA * RB * torch.linalg.norm(B2) / torch.linalg.norm(RB)).requires_grad_()
torch.manual_seed(SEED)
p = [(0.4 * torch.randn(*s, dtype=torch.complex128)).requires_grad_() for s in ((J, D, din), (J, D, din), (P, D), (P, D))]
with torch.no_grad(): X = target(A, B)
r = refit(p, X, 3000)
optA = torch.optim.Adam([A, B], lr=0.01)
best = (r, None)
for it in range(OUTER):
    X = target(A, B)
    with torch.no_grad(): Mp = model(p)
    res = torch.linalg.norm(Mp - X) / torch.linalg.norm(X)
    optA.zero_grad(); (-res).backward(); optA.step()
    with torch.no_grad(): X = target(A, B)
    r = refit(p, X.detach(), 300)
    if r > best[0]: best = (r, X.detach().clone())
    if it % 10 == 0 or it == OUTER - 1:
        ptm = float(torch.linalg.eigvalsh(pt(X.detach()))[0])
        print(f"m = {m}, seed {SEED}, step {it}: residual {r:.3e} (PT-min {ptm:+.4f}); best so far {best[0]:.3e}", flush=True)
# re-check the worst target from fresh random starts, with a larger model
Xw = best[1]
def fresh(X, J, P, RS=3):
    out = np.inf
    for _ in range(RS):
        q = [(0.4 * torch.randn(*s, dtype=torch.complex128)).requires_grad_() for s in ((J, D, din), (J, D, din), (P, D), (P, D))]
        out = min(out, refit(q, X, 4000))
    return out
print(f"m = {m}, seed {SEED}: worst target re-fitted from fresh starts: {J} + {P}: {fresh(Xw, J, P):.3e}; "
      f"{2 * J} + {2 * P}: {fresh(Xw, 2 * J, 2 * P):.3e}", flush=True)
np.save(f"adv_lm6_from_d{DELTA}.npy", Xw.numpy())
