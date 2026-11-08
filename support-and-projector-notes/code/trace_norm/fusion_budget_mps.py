"""Section 9, Conjecture FB on matrix-product inputs (transfer matrices; the merged pair is contracted as one 9-dimensional site): minimise [Delta_full - gamma*Delta_merged]/(s1 s2), first two sites merged.
usage: python fusion_budget_mps.py n chi t gamma runs seed"""
import sys, numpy as np, torch
torch.set_default_dtype(torch.float64)
n, chi, t, gam, runs, seed = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
d = 3
def step(L, As, Bs, swap):
    T1 = torch.einsum('acbd,xae->xcbde', L, As)
    T2 = torch.einsum('xcbde,ybg->xycdeg', T1, Bs)
    if swap: T2 = T2 - t * T2.transpose(0, 1)
    T4 = torch.einsum('uvcdeg,ucf->vdegf', T2, As.conj())
    return torch.einsum('vdegf,vdh->efgh', T4, Bs.conj())
def boundary(l, m):
    return torch.einsum('ia,jc,jb,id->acbd', l, l.conj(), m, m.conj())
def close(L, r, rp):
    return torch.einsum('efgh,e,f,g,h->', L, r, r.conj(), rp, rp.conj())
def merged(Ak, Bk):
    return torch.einsum('xab,ybc->xyac', Ak[0], Ak[1]).reshape(d * d, chi, chi), torch.einsum('xab,ybc->xyac', Bk[0], Bk[1]).reshape(d * d, chi, chi)
def gram(A, l, r):
    L = torch.einsum('ia,jc->ijac', l.conj(), l)
    for s in range(n): L = torch.einsum('ijac,xae,xcf->ijef', L, A[s].conj(), A[s])
    return torch.einsum('ijef,e,f->ij', L, r.conj(), r)
def unpack(ps):
    c = [torch.complex(ps[2 * i], ps[2 * i + 1]) for i in range(len(ps) // 2)]
    return c[:n], c[n:2 * n], c[2 * n], c[2 * n + 1], c[2 * n + 2], c[2 * n + 3]
def evaluate(ps):
    A, B, l, m, r, rp = unpack(ps)
    L0 = boundary(l, m)
    Lf = L0; Ln = L0
    for s in range(n):
        Lf = step(Lf, A[s], B[s], True); Ln = step(Ln, A[s], B[s], False)
    q = close(Lf, r, rp).real; n2 = close(Ln, r, rp).real
    Am, Bm = merged(A[:2], B[:2])
    Lm = step(L0, Am, Bm, True)
    for s in range(2, n): Lm = step(Lm, A[s], B[s], True)
    qm = close(Lm, r, rp).real
    GA, GB = gram(A, l, r), gram(B, m, rp)
    s1s2 = torch.sqrt(torch.clamp(((GA[0,0]*GA[1,1]-GA[0,1]*GA[1,0])*(GB[0,0]*GB[1,1]-GB[0,1]*GB[1,0])).real, min=1e-30))
    Lval = n2 - t * (n2 + 2 * s1s2)
    Df = q - (1 - t) ** (n - 1) * Lval; Dm = qm - (1 - t) ** (n - 2) * Lval
    return Df, Dm, s1s2, n2
rng = np.random.default_rng(seed)
def init():
    ps = []
    for s in range(2 * n):
        for _ in range(2): ps.append(torch.tensor(rng.standard_normal((d, chi, chi)) / np.sqrt(chi), requires_grad=True))
    for shape in ((2, chi), (2, chi), (chi,), (chi,)):
        for _ in range(2): ps.append(torch.tensor(rng.standard_normal(shape), requires_grad=True))
    return ps
print(f"FB on MPS: n={n} chi={chi} t={t} gamma={gam}", flush=True)
best = 1e9
for rr in range(runs):
    ps = init()
    opt = torch.optim.LBFGS(ps, lr=1, max_iter=800, tolerance_grad=1e-12, tolerance_change=1e-15, history_size=50, line_search_fn='strong_wolfe')
    def cl():
        opt.zero_grad(); Df, Dm, ss, n2 = evaluate(ps); v = (Df - gam * Dm) / ss; v.backward(); return v
    for _ in range(3):
        try: opt.step(cl)
        except Exception as e: print('  fail', e); break
    with torch.no_grad(): Df, Dm, ss, n2 = evaluate(ps)
    v = ((Df - gam * Dm) / n2).item()
    if not np.isfinite(v): print(f"run {rr}: non-finite"); continue
    best = min(best, v)
    print(f"run {rr}: (Df - g Dm)/|C|^2 {v:.3e}   Df {Df.item()/n2.item():.3e}  Dm {Dm.item()/n2.item():.3e}  s1s2/|C|^2 {(ss/n2).item():.3f}", flush=True)
print(f"smallest {best:.3e}")
